#!/usr/bin/env python3
"""guidescan_seed.py -- seed-anchored gRNA binding site finder (pipeline stage 1).

Most-likely gRNA binding site per BED interval, found in two stages instead
of one full-peak GuideScan2 enumerate:

  Stage A (grep, no external binary): extract each peak padded with enough
    flanking genome to hold a full spacer+PAM window, then search that padded
    sequence and its reverse complement for occurrences of a required --seed
    (a contiguous substring of --grna + --pam, e.g. the PAM-proximal 8-12 nt
    of the spacer plus the PAM itself, such as ACCAGGACCNGG). 'N' is a
    wildcard in the seed (matches any base) so a seed can straddle or lie
    entirely within the PAM. Matching is otherwise exact -- no mismatches
    within the seed itself. --seed must occur at least once within --grna +
    --pam; if it occurs more than once (e.g. a short/degenerate seed against a
    repetitive spacer), the PAM-proximal occurrence is used, since that's the
    region Cas9 target recognition is most sensitive to. Each hit is
    back-projected to a small candidate window (spacer+PAM length, expanded by
    rna_bulges+dna_bulges nt on each side for bulge headroom).

  Stage B (direct scoring, no external binary): the seed anchor already tells
    us exactly where the candidate spacer+PAM alignment sits, so there is
    nothing left to SEARCH for -- only to SCORE. A small bounded DP (allowing
    up to --rna-bulges/--dna-bulges unpaired bases near that anchor) finds the
    best-scoring placement directly, ranked by: highest bulge-free Doench-2016
    CFD -> fewest bulges -> closest to peak center. This is why it's cheap:
    computing a score for a KNOWN position is O(guide length); GuideScan2's
    enumerate was expensive because it re-searches for the position under a
    wide mismatch/bulge budget, which the seed already made unnecessary.

Rationale: peaks with no exact seed match anywhere are skipped (fast, and
explicit in the output via the seed_hits column); peaks that do have a seed
match get scored directly, in a single tiny DP, with no external tool and no
search over the peak at all.

The Doench-2016 CFD tables come from the sibling doench_cfd.py, shared with
guidescan.py (stage 2) so both stages score identically. No GuideScan2 binary
is needed.

Output TSV columns are chrom/start/end/strand, effective_mismatches,
mismatches, rna_bulges, dna_bulges, cfd, match_sequence, seed_hits (number of
exact-seed anchors found for that peak; 0 means the seed was not found
anywhere in the peak on either strand), then peak_chrom/peak_start/peak_end.

If --raw-hits-outfile is given, also writes the unscored Stage A table: one
row per exact-seed anchor found (multiple rows per peak when more than one
was found), with the raw ungapped gRNA-binding-site and PAM-binding-site
sequences and strand -- exactly what the grep step located, before any DP
scoring, bulge resolution, or per-peak ranking.

Progress: chunks complete out of order (dispatched to a worker pool via
as_completed, not map), so a "progress: X/N chunks (P%), R chunks/s, elapsed
E, ETA T" line is written to stderr roughly every 1% of chunks, using the
observed average throughput so far for the ETA.

Dependencies: stdlib Python and `samtools faidx` only -- no GuideScan2 binary, no
gcc/GLIBCXX workaround. Both come from becaskas.sif via GUIDESCAN_SEED.

Peaks this stage cannot anchor (seed_hits == 0) feed stage 2, GUIDESCAN_ENUMERATE,
which runs the real `enumerate` via the sibling guidescan.py over those and nothing
else. See the README for why the stages are split.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

from doench_cfd import MM_SCORES, PAM_SCORES

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_bed_intervals(bed_path):
    """Return de-duplicated (chrom, start, end) tuples in first-seen order."""
    seen = set()
    intervals = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if len(f) < 3:
                continue
            chrom, start, end = f[0], int(f[1]), int(f[2])
            key = (chrom, start, end)
            if key not in seen:
                seen.add(key)
                intervals.append(key)
    if not intervals:
        sys.exit("ERROR: no usable intervals found in %s" % bed_path)
    return intervals


def _fmt_duration(seconds):
    """H:MM:SS, or 'unknown' for a non-finite (e.g. zero-throughput) ETA."""
    if seconds == float("inf") or seconds != seconds:
        return "unknown"
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    mi, s = divmod(rem, 60)
    return "%d:%02d:%02d" % (h, mi, s)


def parse_args():
    p = argparse.ArgumentParser(
        description="Most-likely gRNA binding site per BED interval via exact-seed grep "
                     "+ direct DP scoring + CFD.")
    p.add_argument("--bed", required=True, help="Input BED (uses cols 1-3; extras ignored).")
    p.add_argument("--genome", required=True, help="Indexed FASTA (.fai present; UCSC chr naming).")
    p.add_argument("--grna", required=True, help="Spacer sequence, e.g. CCCTACCTGTCACCAGGACC.")
    p.add_argument("--pam", default="NGG", help="PAM pattern (default NGG; N=wildcard).")
    p.add_argument("--seed", required=True,
                   help="Seed substring of --grna + --pam, grepped verbatim (0 mismatches "
                        "except at 'N' positions, which are wildcards matching any base -- "
                        "use this to include the PAM, e.g. ACCAGGACCNGG). Must occur at least "
                        "once in --grna + --pam; if it occurs more than once, the PAM-proximal "
                        "occurrence is used.")
    p.add_argument("--mismatches", type=int, default=10,
                   help="Discard a scored site with more than this many mismatches "
                        "(default 10). Applied AFTER scoring, as a filter -- not a search "
                        "budget, since Stage B no longer searches.")
    p.add_argument("--rna-bulges", type=int, default=2, dest="rna_bulges",
                   help="Max RNA bulges to consider near the anchor (default 2).")
    p.add_argument("--dna-bulges", type=int, default=2, dest="dna_bulges",
                   help="Max DNA bulges to consider near the anchor (default 2).")
    p.add_argument("--samtools", default="samtools", help="Path to samtools (default: on PATH).")
    p.add_argument("--chunk-size", type=int, default=15, dest="chunk_size",
                   help="Peaks processed together as one padded-extraction unit. Default 15.")
    p.add_argument("--summit-width", type=int, default=0, dest="summit_width",
                   help="If >0, restrict the Stage A seed search to the center N bp of each "
                        "peak instead of its full BED span (plus the usual bulge-headroom "
                        "margin around that narrowed region). Output still reports each peak's "
                        "original full BED bounds. Default 0 (search the full peak).")
    p.add_argument("--workers", type=int, default=(os.cpu_count() or 1),
                   help="Chunks to process in parallel (default: all CPUs).")
    p.add_argument("--outfile", required=True, help="Output TSV path.")
    p.add_argument("--raw-hits-outfile", default=None, dest="raw_hits_outfile",
                   help="If given, also write the raw Stage A table: one row per exact-seed "
                        "anchor (multiple rows per peak when more than one was found), with "
                        "the ungapped gRNA-binding-site and PAM-binding-site sequences and "
                        "strand, before any DP scoring or ranking.")
    p.add_argument("--workdir", default=None,
                   help="Parent dir for per-chunk temp dirs (default: system temp; "
                        "point at $L_SCRATCH/$SCRATCH for very large runs).")
    return p.parse_args()


def summit_region(key, width):
    """Narrow peak `key` (chrom, s, e) to its center `width` bp; returns `key`
    unchanged if width <= 0. Odd `width` puts the extra base downstream.
    """
    if width <= 0:
        return key
    chrom, s, e = key
    center = (s + e) // 2
    ss = center - width // 2
    return (chrom, max(0, ss), ss + width)


def extract_padded_sequences(intervals, margin, genome, samtools, workdir, summit_width=0):
    """Extract each peak's search region -- its center `summit_width` bp if
    summit_width > 0, else its full (chrom, start, end) span -- padded by
    `margin` nt each side.

    Returns {interval: (chrom, pad_start, pad_end, seq)}, keyed by the
    ORIGINAL (unnarrowed) interval; pad_start/pad_end are 0-based half-open
    genomic coords (start clamped at 0), seq is uppercase.
    """
    # Ordered (not dict-keyed) on purpose: narrowing every peak to the same fixed
    # summit_width means two DIFFERENT peaks can share an identical midpoint and
    # therefore an identical padded region string. samtools faidx -r emits one
    # record per region-file line in order (duplicates and all), so we track
    # peaks by position instead of by region string to avoid silently colliding
    # two peaks onto one dict entry.
    regfile = os.path.join(workdir, "padded_regions.txt")
    order = []
    with open(regfile, "w") as fh:
        for key in intervals:
            chrom, s, e = summit_region(key, summit_width)
            ps, pe = max(0, s - margin), e + margin
            order.append((key, chrom, ps, pe))
            fh.write("%s:%d-%d\n" % (chrom, ps + 1, pe))

    raw = subprocess.run([samtools, "faidx", genome, "-r", regfile],
                         check=True, capture_output=True, text=True).stdout

    seqs = []
    name_seen = False
    buf = []
    for line in raw.splitlines():
        if line.startswith(">"):
            if name_seen:
                seqs.append("".join(buf))     # flush even if empty, to stay aligned with `order`
            name_seen = True
            buf = []
        else:
            buf.append(line.strip())
    if name_seen:
        seqs.append("".join(buf))

    if len(seqs) != len(order):
        sys.exit("ERROR: samtools faidx returned %d record(s) for %d requested region(s)"
                 % (len(seqs), len(order)))

    out = {}
    for (key, chrom, ps, pe), seq in zip(order, seqs):
        if not seq:
            sys.stderr.write("WARNING: no sequence returned for %s:%d-%d (skipped)\n"
                             % (chrom, ps + 1, pe))
        out[key] = (chrom, ps, pe, seq.upper())
    return out


_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT",
}


def _iupac_eq(a, b):
    """True if IUPAC codes `a`/`b` are compatible -- their possible-base sets
    intersect (e.g. 'N' matches anything; 'W' [A or T] matches 'A' or 'T' but
    not 'C'/'G'). Unrecognized characters are treated as literal (compared
    for exact equality).
    """
    return bool(set(_IUPAC.get(a, a)) & set(_IUPAC.get(b, b)))


def find_pattern(seq, pattern):
    """Yield every start index where `pattern` matches `seq`, treating IUPAC
    ambiguity codes (N, W, R, Y, ...) in either string as their possible-base
    set rather than a literal character. O(len(seq) * len(pattern)) -- fine
    here since it only ever runs over a single padded peak / a spacer+PAM,
    not the whole genome.
    """
    n, m = len(seq), len(pattern)
    for i in range(n - m + 1):
        if all(_iupac_eq(pattern[k], seq[i + k]) for k in range(m)):
            yield i


def find_seed_anchors(target, seed, seed_offset):
    """Exact-seed hits in `target` ('N' in the seed is a wildcard), back-projected
    to the local start of the ideal (0-bulge) spacer alignment. `target` is
    whatever strand's sequence this anchor should be scored against (the
    caller passes `seq` for '+' hits, its reverse complement for '-' hits).
    Returns a sorted, deduped list of anchor positions (may be out of bounds;
    score_local_alignment checks that).
    """
    return sorted({i - seed_offset for i in find_pattern(target, seed)})


def _better(a, b):
    """True if DP value `a` = (edits, score, mismatches) beats `b` (or `b` is
    None): fewer edits wins first; CFD-product is only a tiebreak among
    equally-parsimonious alignments -- bulges must never be a free escape
    from a mismatch penalty, or the DP would invent bulges purely to inflate
    CFD instead of reporting what's actually the most parsimonious alignment.
    """
    if b is None:
        return True
    if a[0] != b[0]:
        return a[0] < b[0]
    return a[1] > b[1]


def score_local_alignment(target, anchor, guide, pam, max_rna, max_dna):
    """Minimum-edit alignment of `guide`+`pam` against `target`, anchored at
    `anchor` (the seed-derived, ideal 0-bulge start), allowing up to `max_rna`
    RNA bulges (unpaired guide bases) and `max_dna` DNA bulges (unpaired
    target bases) within the spacer; CFD is the tiebreak among equally-
    parsimonious alignments, not the primary objective. This SCORES a known
    position -- it does not search for one -- so it's a small bounded DP
    (O(len(guide) * max_rna * max_dna), a few hundred states), not a search.

    Returns a dict (mm, rna, dna, cfd, start, end) for the best placement, or
    None if `anchor` doesn't leave room for a full site.
    """
    n, plen, tlen = len(guide), len(pam), len(target)
    dp = {(0, 0): (0, 1.0, 0)}   # (r, d) -> (edits, score, mismatches) at guide position i
    for i in range(n):
        # DNA bulges don't advance i, so saturate same-i transitions first
        # (bounded: d only increases, capped at max_dna).
        for _ in range(max_dna):
            grown = {}
            for (r, d), (ed, score, mm) in dp.items():
                if d >= max_dna:
                    continue
                tpos = anchor + (i + d - r)
                if 0 <= tpos < tlen:
                    key, cand = (r, d + 1), (ed + 1, score, mm)
                    if _better(cand, dp.get(key)) and _better(cand, grown.get(key)):
                        grown[key] = cand
            if not grown:
                break
            for k, v in grown.items():
                if _better(v, dp.get(k)):
                    dp[k] = v

        nxt = {}
        gb = guide[i]
        r_key = "U" if gb == "T" else gb
        for (r, d), (ed, score, mm) in dp.items():
            tpos = anchor + (i + d - r)
            if 0 <= tpos < tlen:
                tb = target[tpos]
                if tb == gb:
                    cand = (ed, score, mm)
                else:
                    dbase = tb.translate(COMP)
                    factor = MM_SCORES.get("r%s:d%s,%d" % (r_key, dbase, i + 1), 1.0)
                    cand = (ed + 1, score * factor, mm + 1)
                key = (r, d)
                if _better(cand, nxt.get(key)):
                    nxt[key] = cand
            if r < max_rna:
                key, cand = (r + 1, d), (ed + 1, score, mm)
                if _better(cand, nxt.get(key)):
                    nxt[key] = cand
        dp = nxt
        if not dp:
            return None

    best = None
    for (r, d), (ed, score, mm) in dp.items():
        end = anchor + (n + d - r)
        if end < 0 or end + plen > tlen:
            continue
        pam_factor = PAM_SCORES.get(target[end:end + plen][1:3], 0.0)
        cfd = score * pam_factor
        rank_key = (ed, -cfd)
        if best is None or rank_key < best[0]:
            best = (rank_key, dict(mm=mm, rna=r, dna=d, cfd=cfd, start=anchor, end=end))
    return best[1] if best else None


def raw_hit_row(peak_key, strand, anchor, target, pad_chrom, pad_start, n_seq, guide_len, pam_len):
    """The ungapped gRNA-binding-site and PAM-binding-site sequences at one
    Stage A anchor, exactly as the seed located them -- no DP, no bulge
    resolution, no ranking. Returns None if the anchor doesn't leave room for
    a full site (mirrors score_local_alignment's own bounds check).
    """
    end_spacer = anchor + guide_len
    end_pam = end_spacer + pam_len
    if anchor < 0 or end_pam > len(target):
        return None
    if strand == "+":
        site_start, site_end = pad_start + anchor, pad_start + end_pam
    else:
        site_start, site_end = pad_start + (n_seq - end_pam), pad_start + (n_seq - anchor)
    _pchrom, ps, pe = peak_key
    return dict(peak_chrom=pad_chrom, peak_start=ps, peak_end=pe, strand=strand,
                grna_binding_site=target[anchor:end_spacer],
                pam_binding_site=target[end_spacer:end_pam],
                site_chrom=pad_chrom, site_start=site_start, site_end=site_end)


def process_chunk(job):
    """Stage A (seed anchors) + Stage B (direct DP scoring) for one chunk of
    peaks -- no external binary, no search: the seed already tells us where
    to look, so Stage B only ever scores that one anchored placement per
    strand. Returns (best, seed_hits, raw_rows): `best` maps peak -> ranked
    hit dict (peaks with a scorable anchor -- see `mismatches` filter below);
    `seed_hits` maps every peak in the chunk -> count of exact-seed anchors
    found (0 if none, on either strand); `raw_rows` is one row per anchor
    (multiple rows per peak when Stage A found more than one), from
    `raw_hit_row` above -- the raw grep-step table, before any scoring.

    `job` is a plain tuple so it pickles cleanly to ProcessPoolExecutor workers.
    """
    (chunk, genome, samtools, guide, pam, seed, seed_offset,
     mismatches, rna_b, dna_b, summit_width, parent_tmp) = job
    wd = tempfile.mkdtemp(prefix="grna_seed_chunk_", dir=parent_tmp)
    try:
        spacer_len, pam_len = len(guide), len(pam)
        margin = spacer_len + pam_len + rna_b + dna_b
        padded = extract_padded_sequences(chunk, margin, genome, samtools, wd, summit_width)

        best, seed_hits, raw_rows = {}, {}, []
        for peak_key in chunk:
            pad_chrom, pad_start, _pad_end, seq = padded[peak_key]
            n_seq = len(seq)
            rc = seq.translate(COMP)[::-1]

            candidates = [("+", a) for a in find_seed_anchors(seq, seed, seed_offset)]
            candidates += [("-", a) for a in find_seed_anchors(rc, seed, seed_offset)]
            seed_hits[peak_key] = len(candidates)

            _pchrom, ps, pe = peak_key
            peak_center = ps + (pe - ps) / 2.0
            best_rec = None
            for strand, anchor in candidates:
                target = seq if strand == "+" else rc
                row = raw_hit_row(peak_key, strand, anchor, target, pad_chrom, pad_start,
                                   n_seq, spacer_len, pam_len)
                if row is not None:
                    raw_rows.append(row)
                r = score_local_alignment(target, anchor, guide, pam, rna_b, dna_b)
                if r is None or r["mm"] > mismatches:
                    continue
                if strand == "+":
                    site_start, site_end = pad_start + r["start"], pad_start + r["end"]
                else:
                    site_start = pad_start + (n_seq - r["end"])
                    site_end = pad_start + (n_seq - r["start"])
                cfd = r["cfd"]
                bulges = r["rna"] + r["dna"]
                site_center = site_start + (site_end - site_start) / 2.0
                cdist = abs(site_center - peak_center)
                cfd_rank = cfd if cfd == cfd else -1.0
                key = (-cfd_rank, bulges, cdist)
                if best_rec is None or key < best_rec["_key"]:
                    best_rec = dict(chrom=pad_chrom, start=site_start, end=site_end,
                                     strand=strand, mm=r["mm"], rna=r["rna"], dna=r["dna"],
                                     eff_mm=r["mm"] + bulges, cfd=cfd,
                                     seq=target[r["start"]:r["end"] + pam_len],
                                     pchrom=pad_chrom, pstart=ps, pend=pe, _key=key)
            if best_rec is not None:
                best[peak_key] = best_rec
        return best, seed_hits, raw_rows
    finally:
        shutil.rmtree(wd, ignore_errors=True)


def main():
    args = parse_args()
    guide = args.grna.upper()
    pam = args.pam.upper()
    seed = args.seed.upper()

    if shutil.which(args.samtools) is None and not os.path.exists(args.samtools):
        sys.exit("ERROR: samtools not found: %s" % args.samtools)
    if not os.path.exists(args.genome + ".fai"):
        sys.exit("ERROR: FASTA index missing: %s.fai (run `samtools faidx`)" % args.genome)
    if args.workdir:
        os.makedirs(args.workdir, exist_ok=True)

    target = guide + pam    # seed may lie within the spacer, straddle into the PAM, or both
    occurrences = list(find_pattern(target, seed))
    if not occurrences:
        sys.exit("ERROR: --seed %r not found in --grna+--pam %r ('N' in --seed is a wildcard)"
                 % (seed, target))
    seed_offset = max(occurrences)     # ambiguous seed -> always prefer the PAM-proximal match
    if len(occurrences) > 1:
        sys.stderr.write("WARNING: --seed %r matches --grna+--pam %r at %d positions %r; "
                         "using the PAM-proximal one (offset %d)\n"
                         % (seed, target, len(occurrences), occurrences, seed_offset))
    intervals = read_bed_intervals(args.bed)
    cs = max(1, args.chunk_size)
    chunks = [intervals[i:i + cs] for i in range(0, len(intervals), cs)]
    jobs = [(chunk, args.genome, args.samtools,
             guide, pam, seed, seed_offset, args.mismatches, args.rna_bulges, args.dna_bulges,
             args.summit_width, args.workdir)
            for chunk in chunks]
    workers = max(1, min(args.workers, len(jobs)))
    sys.stderr.write("%d unique interval(s) -> %d chunk(s) of %d, %d parallel worker(s), seed=%r "
                     "(offset %d), summit_width=%d\n"
                     % (len(intervals), len(chunks), cs, workers, seed, seed_offset, args.summit_width))

    total_jobs = len(jobs)
    report_every = max(1, total_jobs // 100)   # ~100 progress lines across the whole run
    t_start = time.time()

    def report_progress(done):
        if done != total_jobs and done % report_every != 0:
            return
        elapsed = time.time() - t_start
        rate = done / elapsed if elapsed > 0 else 0.0
        eta = (total_jobs - done) / rate if rate > 0 else float("inf")
        sys.stderr.write("progress: %d/%d chunks (%.1f%%), %.2f chunks/s, elapsed %s, ETA %s\n"
                         % (done, total_jobs, 100.0 * done / total_jobs, rate,
                            _fmt_duration(elapsed), _fmt_duration(eta)))

    results = []
    if workers == 1 or len(jobs) == 1:
        for i, j in enumerate(jobs, 1):
            results.append(process_chunk(j))
            report_progress(i)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = [ex.submit(process_chunk, j) for j in jobs]
            for i, fut in enumerate(as_completed(futures), 1):
                results.append(fut.result())
                report_progress(i)

    best, seed_hits, raw_rows = {}, {}, []
    for b, sh, rr in results:
        best.update(b)
        seed_hits.update(sh)
        raw_rows.extend(rr)

    if args.raw_hits_outfile:
        raw_cols = ["peak_chrom", "peak_start", "peak_end", "strand",
                    "grna_binding_site", "pam_binding_site",
                    "site_chrom", "site_start", "site_end"]
        raw_rows.sort(key=lambda r: (r["peak_chrom"], r["peak_start"], r["strand"], r["site_start"]))
        with open(args.raw_hits_outfile, "w") as out:
            out.write("\t".join(raw_cols) + "\n")
            for r in raw_rows:
                out.write("\t".join(str(r[c]) for c in raw_cols) + "\n")
        sys.stderr.write("Wrote %s: %d raw seed anchor(s) across %d peak(s)\n"
                         % (args.raw_hits_outfile, len(raw_rows),
                            sum(1 for v in seed_hits.values() if v > 0)))

    cols = ["chrom", "start", "end", "strand", "effective_mismatches", "mismatches",
            "rna_bulges", "dna_bulges", "cfd", "match_sequence", "seed_hits",
            "peak_chrom", "peak_start", "peak_end"]
    n_hit = 0
    with open(args.outfile, "w") as out:
        out.write("\t".join(cols) + "\n")
        for peak_key in intervals:      # original BED order
            chrom, s, e = peak_key
            nseed = seed_hits.get(peak_key, 0)
            if peak_key in best:
                r = best[peak_key]
                cfd_str = "NA" if r["cfd"] != r["cfd"] else "%.6g" % r["cfd"]
                out.write("\t".join(str(x) for x in [
                    r["chrom"], r["start"], r["end"], r["strand"], r["eff_mm"], r["mm"],
                    r["rna"], r["dna"], cfd_str, r["seq"], nseed,
                    r["pchrom"], r["pstart"], r["pend"]]) + "\n")
                n_hit += 1
            else:
                out.write("\t".join(str(x) for x in [chrom, ".", ".", ".", "NA", "NA", "NA", "NA",
                                     "NA", "NA", nseed, chrom, s, e]) + "\n")
    n_no_seed = sum(1 for k in intervals if seed_hits.get(k, 0) == 0)
    sys.stderr.write("Wrote %s: %d peak(s), %d with a scored site, %d with 0 seed hits\n"
                     % (args.outfile, len(intervals), n_hit, n_no_seed))



if __name__ == "__main__":
    main()
