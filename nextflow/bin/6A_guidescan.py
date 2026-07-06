#!/usr/bin/env python3
"""grna_binding_sites.py -- most-likely gRNA binding site per BED interval, via GuideScan2.

Given a BED file and a gRNA (spacer + PAM), this runs GuideScan2 `enumerate` over a
peaks-only FM-index built from the BED intervals, then reports, per interval, the single
most likely protospacer ranked by:  highest bulge-free CFD -> fewest bulges ->
closest to peak center.  CFD (Doench 2016) is a product of per-position mismatch
penalties, so it already encodes every mismatch; mismatches are therefore NOT
scored a second time in the ranking (that would double-count them).  Bulges, which
CFD cannot score, are the tiebreak.  Output is a headered TSV whose first three
columns form a valid BED of the predicted binding site (source-peak coordinates
kept as extra columns for provenance).

CFD is the standard Doench et al. 2016 score, computed in-script from embedded tables
(defined at the bottom of this file, copied verbatim from GuideScan2's doench.hpp), so
the script needs only a stock GuideScan2 binary + samtools at runtime -- GuideScan2 is
never modified.  It is computed "bulge-free": an RNA bulge (unpaired guide base, marked
'.' by GuideScan's complete mode) contributes a neutral factor (that guide position has
no DNA partner) and a DNA bulge (an inserted DNA base) is removed before scoring, so a
bulge hit still gets a real CFD reflecting its mismatch pattern.  Guide positions are
intrinsic (1..n) and are NOT renumbered around bulges.

Coordinates: BED is 0-based half-open; samtools faidx is 1-based inclusive; GuideScan2
reports the match's 1-based lower coordinate within the record.  Site 0-based genomic
start = peak_bedstart + (match_position - 1); end = start + len(spacer) + len(PAM).

Parallelism: peaks are split into chunks of --chunk-size (default 15), and each chunk is
processed as an independent unit (its own tiny index + enumerate) across a pool of
--workers processes. GuideScan2's enumerate is single-threaded per guide, so chunking
across peaks -- not adding threads -- is what actually scales this to many peaks.

Dependencies (this script is a self-contained single file, but it orchestrates two
external executables and needs an indexed genome -- it is NOT dependency-free at runtime):

  * Python           >= 3.6   (stdlib only: argparse, os, shutil, subprocess, sys,
                                tempfile; no third-party packages). On Sherlock:
                                `ml python/3.12.1`.
  * samtools                  used for `faidx` sequence extraction; needs the genome
                                FASTA index (<genome>.fai). Tested with 1.16.1.
                                (htslib links libcurl at runtime -- e.g. curl 8.x.)
  * GuideScan2 (guidescan)    the `index` + `enumerate` subcommands. Stock, UNMODIFIED
                                binary. Tested with GuideScan2 2.0.0
                                (github.com/pritykinlab/guidescan-cli). If built with
                                GCC >= 13, add `#include <cstdint>` to
                                include/genomics/structures.hpp before compiling.
  * Genome FASTA              indexed (`samtools faidx genome.fasta` -> .fai present),
                                UCSC-style chromosome names (chr1, chr2, ...). Tested
                                with hg38.

Sherlock environment setup (load before running):
  ml python/3.12.1 biology samtools/1.16.1 system curl/8.17.0
  # curl module is required -- samtools/htslib and the guidescan binary link libcurl.

The embedded CFD tables (bottom of file) are the Doench et al. 2016 sgRNA mismatch/PAM
scores, copied verbatim from GuideScan2's include/genomics/doench.hpp.

Paths to `samtools` and `guidescan` are passed via --samtools / --guidescan (default:
found on PATH); the script fails fast with a clear error if either is missing.
"""
import argparse
import itertools
import os
import shutil
import subprocess
import sys
import tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_args():
    p = argparse.ArgumentParser(
        description="Most-likely gRNA binding site per BED interval via GuideScan2 + CFD.")
    p.add_argument("--bed", required=True, help="Input BED (uses cols 1-3; extras ignored).")
    p.add_argument("--genome", required=True, help="Indexed FASTA (.fai present; UCSC chr naming).")
    p.add_argument("--grna", required=True, help="Spacer sequence, e.g. CCCTACCTGTCACCAGGACC.")
    p.add_argument("--pam", default="NGG", help="PAM pattern (default NGG; N=wildcard).")
    p.add_argument("--mismatches", type=int, default=10, help="Max mismatches in the spacer (default 10).")
    p.add_argument("--rna-bulges", type=int, default=2, dest="rna_bulges", help="Max RNA bulges (default 2).")
    p.add_argument("--dna-bulges", type=int, default=2, dest="dna_bulges", help="Max DNA bulges (default 2).")
    p.add_argument("--guidescan", default="guidescan", help="Path to guidescan binary (default: on PATH).")
    p.add_argument("--samtools", default="samtools", help="Path to samtools (default: on PATH).")
    p.add_argument("--chunk-size", type=int, default=15, dest="chunk_size",
                   help="Peaks processed together as one unit (own index+enumerate). Default 15.")
    p.add_argument("--workers", type=int, default=(os.cpu_count() or 1),
                   help="Chunks to process in parallel (default: all CPUs).")
    p.add_argument("--threads", type=int, default=1,
                   help="guidescan -n per chunk (single guide, so leave at 1).")
    p.add_argument("--outfile", required=True, help="Output TSV path.")
    p.add_argument("--workdir", default=None,
                   help="Parent dir for per-chunk temp dirs (default: system temp; "
                        "point at $L_SCRATCH/$SCRATCH for very large runs).")
    return p.parse_args()


def run(cmd):
    """Run a command quietly; on failure raise with captured output (tails)."""
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("command failed: %s\n%s\n%s"
                           % (" ".join(str(c) for c in cmd), r.stdout[-2000:], r.stderr[-2000:]))


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


def peak_index(name):
    """Global 1-based index encoded in a peak record name (peak000123_chr_..)."""
    return int(name.split("_", 1)[0][4:])


def build_peaks_fasta(intervals, genome, samtools, workdir, start_index=1):
    """Extract each interval and write a FASTA whose record names encode coordinates.

    Record name = peak{globalindex}_{chrom}_{bedstart}_{bedend}; indices start at
    `start_index` so names stay globally unique/ordered across parallel chunks.
    Returns (peaks_fa, {name: (chrom, s, e)}).
    """
    regfile = os.path.join(workdir, "regions.txt")
    # faidx region = 1-based inclusive: (bedstart+1)..bedend
    reg_to_meta = {}
    with open(regfile, "w") as fh:
        for i, (chrom, s, e) in enumerate(intervals, start_index):
            reg = "%s:%d-%d" % (chrom, s + 1, e)
            reg_to_meta[reg] = (i, chrom, s, e)
            fh.write(reg + "\n")

    raw = subprocess.run([samtools, "faidx", genome, "-r", regfile],
                         check=True, capture_output=True, text=True).stdout

    # Parse faidx multi-FASTA; headers are the region strings we requested.
    seqs = {}
    name = None
    buf = []
    for line in raw.splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)

    peaks_fa = os.path.join(workdir, "peaks.fa")
    name_meta = {}
    n_written = 0
    with open(peaks_fa, "w") as fh:
        for reg, (idx, chrom, s, e) in reg_to_meta.items():
            if reg not in seqs or not seqs[reg]:
                sys.stderr.write("WARNING: no sequence returned for %s (skipped)\n" % reg)
                continue
            rec = "peak%06d_%s_%d_%d" % (idx, chrom, s, e)
            name_meta[rec] = (chrom, s, e)
            fh.write(">%s\n%s\n" % (rec, seqs[reg]))
            n_written += 1
    if n_written == 0:
        sys.exit("ERROR: no peak sequences extracted; check --genome / .fai / chrom names.")
    return peaks_fa, name_meta


def write_kmers(grna, pam, workdir):
    kmers = os.path.join(workdir, "kmers.csv")
    with open(kmers, "w") as fh:
        fh.write("id,sequence,pam,chromosome,position,sense\n")
        fh.write("query,%s,%s,chrNA,1,+\n" % (grna.upper(), pam.upper()))
    return kmers


def cfd_score(guide, match_seq, pam_len, n_dna_bulges):
    """Bulge-free Doench-2016 CFD for a (possibly gapped) hit.

    `match_seq` is GuideScan2 complete-mode output in guide orientation: the spacer
    region followed by the PAM, where matched bases are UPPERCASE, mismatched bases
    are lowercase, RNA bulges (unpaired guide base) are '.', and DNA bulges (an
    inserted DNA base that consumes no guide position) are extra lowercase bases in
    the spacer.  We score CFD "as if the bulges were not there":

      * RNA bulge ('.') -> that guide position has no DNA partner; neutral factor 1.0.
      * DNA bulge        -> an inserted DNA base; removed before scoring.
      * mismatch         -> the standard per-position CFD factor at its guide position.

    Guide positions are intrinsic (1..n, PAM-distal..proximal, matching GuideScan's own
    calculate_cfd) and are NOT renumbered around bulges.  Perfect NGG match -> 1.0.

    Ambiguity note: GuideScan collapses mismatch and DNA-bulge bases to lowercase, so
    for a hit with >=1 DNA bulge we cannot tell which lowercase base is the insertion.
    We return the MAX CFD over all consistent placements (best-case complementarity);
    bulges are penalised separately as the ranking tiebreak.  With 0 DNA bulges the
    parse is exact.
    """
    n = len(guide)
    if pam_len < 3:
        return float("nan")
    spacer = match_seq[:-pam_len]
    pam = match_seq[-pam_len:].upper()
    pam_factor = PAM_SCORES.get(pam[1:3], 0.0)

    # Extra chars beyond the n guide positions are the DNA-bulge insertions.
    n_dna_bulges = max(0, len(spacer) - n) if n_dna_bulges else 0
    lower_idx = [i for i, ch in enumerate(spacer) if ch != "." and ch.islower()]

    def walk(removed):
        """CFD over the spacer with the chars in `removed` (DNA bulges) excised."""
        score, g = 1.0, 0
        for i, ch in enumerate(spacer):
            if i in removed:
                continue                       # DNA bulge: not a guide position
            if g >= n:
                break
            if ch == ".":                      # RNA bulge: no DNA partner, factor 1.0
                g += 1
                continue
            gb = guide[g]
            if ch.upper() != gb:
                r = "U" if gb == "T" else gb
                d = ch.upper().translate(COMP)
                score *= MM_SCORES.get("r%s:d%s,%d" % (r, d, g + 1), 1.0)
            g += 1
        return score

    if n_dna_bulges == 0:
        best = walk(set())
    else:
        best = 0.0
        for combo in itertools.combinations(lower_idx, min(n_dna_bulges, len(lower_idx))):
            best = max(best, walk(set(combo)))
    return best * pam_factor


def rank_hits(enum_csv, name_meta, guide, pam):
    """Stream enumerate CSV; keep the best hit per peak.

    Best = (highest bulge-free CFD, fewest bulges, smallest distance to peak center).
    CFD (Doench-2016) is a product of per-position mismatch penalties, so it already
    encodes every mismatch -- mismatches are NOT scored again in the key (no
    double-counting). Bulges, which CFD cannot represent, are the tiebreak.
    eff_mm is still reported for provenance but does not drive the ranking.
    Returns {peak_name: best_record_dict}.
    """
    n, plen = len(guide), len(pam)
    best = {}
    with open(enum_csv) as fh:
        header = fh.readline().rstrip("\n").split(",")
        idx = {c: i for i, c in enumerate(header)}
        ci_chrm, ci_pos, ci_str = idx["match_chrm"], idx["match_position"], idx["match_strand"]
        ci_dist, ci_seq = idx["match_distance"], idx["match_sequence"]
        ci_rna, ci_dna = idx["rna_bulges"], idx["dna_bulges"]
        for line in fh:
            f = line.rstrip("\n").split(",")
            peak = f[ci_chrm]
            if peak not in name_meta:      # skip NA / no-match rows
                continue
            mm = int(f[ci_dist])
            rna = int(f[ci_rna])
            dna = int(f[ci_dna])
            eff_mm = mm + rna + dna         # provenance only; NOT used to rank
            bulges = rna + dna
            pos = int(f[ci_pos])
            mseq = f[ci_seq]
            # Bulge-free CFD: mismatches are fully captured here, so they are not
            # counted independently in the key; bulges break ties (see below).
            score = cfd_score(guide, mseq, plen, dna)
            cfd_rank = score if score == score else -1.0   # NaN (unscoreable) -> worst
            chrom, bstart, bend = name_meta[peak]
            site_start = bstart + (pos - 1)
            # genomic footprint: guide positions - RNA bulges (no ref base)
            #                    + DNA bulges (extra ref base) + PAM. Approximate
            #                    for bulge hits (guidescan's start assumes ungapped width).
            span = n - rna + dna + plen
            site_end = site_start + span
            site_center = site_start + span / 2.0
            peak_center = bstart + (bend - bstart) / 2.0
            cdist = abs(site_center - peak_center)
            key = (-cfd_rank, bulges, cdist)      # highest CFD, then fewest bulges, then closest
            rec = dict(peak=peak, chrom=chrom, start=site_start, end=site_end,
                       strand=f[ci_str], mm=mm, rna=rna, dna=dna, eff_mm=eff_mm,
                       cfd=score, seq=mseq, pchrom=chrom, pstart=bstart, pend=bend, _key=key)
            if peak not in best or key < best[peak]["_key"]:
                best[peak] = rec
    return best


def process_chunk(job):
    """Process one chunk of peaks as an independent unit: own temp dir, own tiny
    GuideScan2 index + enumerate, own ranking. Returns (best, name_meta) for the chunk.

    `job` is a plain tuple so it pickles cleanly to ProcessPoolExecutor workers.
    """
    (chunk, start_index, genome, samtools, guidescan, guide, pam,
     mismatches, rna_b, dna_b, threads, parent_tmp) = job
    wd = tempfile.mkdtemp(prefix="grna_chunk_", dir=parent_tmp)
    try:
        peaks_fa, name_meta = build_peaks_fasta(chunk, genome, samtools, wd, start_index)
        kmers = write_kmers(guide, pam, wd)
        idxp = os.path.join(wd, "idx")
        enum = os.path.join(wd, "enum.csv")
        run([guidescan, "index", "--index", idxp, peaks_fa])
        run([guidescan, "enumerate",
             "-m", mismatches, "--rna-bulges", rna_b, "--dna-bulges", dna_b,
             "--format", "csv", "--mode", "complete", "-n", threads,
             "-f", kmers, "-o", enum, idxp])
        best = rank_hits(enum, name_meta, guide, pam)
        return best, name_meta
    finally:
        shutil.rmtree(wd, ignore_errors=True)


def main():
    args = parse_args()
    guide = args.grna.upper()
    pam = args.pam.upper()

    for tool, name in ((args.guidescan, "guidescan"), (args.samtools, "samtools")):
        if shutil.which(tool) is None and not os.path.exists(tool):
            sys.exit("ERROR: %s not found: %s" % (name, tool))
    if not os.path.exists(args.genome + ".fai"):
        sys.exit("ERROR: FASTA index missing: %s.fai (run `samtools faidx`)" % args.genome)
    if args.workdir:
        os.makedirs(args.workdir, exist_ok=True)

    intervals = read_bed_intervals(args.bed)
    cs = max(1, args.chunk_size)
    chunks = [intervals[i:i + cs] for i in range(0, len(intervals), cs)]
    jobs = [(chunk, ci * cs + 1, args.genome, args.samtools, args.guidescan,
             guide, pam, args.mismatches, args.rna_bulges, args.dna_bulges,
             args.threads, args.workdir)
            for ci, chunk in enumerate(chunks)]
    workers = max(1, min(args.workers, len(jobs)))
    sys.stderr.write("%d unique interval(s) -> %d chunk(s) of %d, %d parallel worker(s)\n"
                     % (len(intervals), len(chunks), cs, workers))

    # Run chunks (serial if a single worker/chunk; otherwise a process pool).
    if workers == 1 or len(jobs) == 1:
        results = [process_chunk(j) for j in jobs]
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(process_chunk, jobs))

    best, name_meta = {}, {}
    for b, nm in results:
        best.update(b)
        name_meta.update(nm)

    cols = ["chrom", "start", "end", "strand", "effective_mismatches", "mismatches",
            "rna_bulges", "dna_bulges", "cfd", "match_sequence",
            "peak_chrom", "peak_start", "peak_end"]
    n_hit = 0
    with open(args.outfile, "w") as out:
        out.write("\t".join(cols) + "\n")
        # one row per input interval, in BED order; NA when no hit within budget
        for peak, (chrom, s, e) in sorted(name_meta.items(), key=lambda kv: peak_index(kv[0])):
            if peak in best:
                r = best[peak]
                cfd_str = "NA" if r["cfd"] != r["cfd"] else "%.6g" % r["cfd"]
                out.write("\t".join(str(x) for x in [
                    r["chrom"], r["start"], r["end"], r["strand"], r["eff_mm"], r["mm"],
                    r["rna"], r["dna"], cfd_str, r["seq"],
                    r["pchrom"], r["pstart"], r["pend"]]) + "\n")
                n_hit += 1
            else:
                out.write("\t".join([chrom, ".", ".", ".", "NA", "NA", "NA", "NA",
                                     "NA", "NA", chrom, str(s), str(e)]) + "\n")
    sys.stderr.write("Wrote %s: %d peak(s), %d with a site <=%d effective mm (%d NA)\n"
                     % (args.outfile, len(name_meta), n_hit, args.mismatches,
                        len(name_meta) - n_hit))


# ---- Embedded Doench et al. 2016 CFD tables (verbatim from GuideScan2 doench.hpp) ----
MM_SCORES = {
    'rU:dT,12': 0.8,
    'rU:dT,13': 0.692307692,
    'rU:dC,5': 0.64,
    'rG:dA,14': 0.26666666699999997,
    'rG:dG,19': 0.448275862,
    'rG:dG,18': 0.47619047600000003,
    'rG:dG,15': 0.272727273,
    'rG:dG,14': 0.428571429,
    'rG:dG,17': 0.235294118,
    'rG:dG,16': 0.0,
    'rC:dC,20': 0.058823529000000006,
    'rG:dT,20': 0.9375,
    'rG:dG,13': 0.42105263200000004,
    'rG:dG,12': 0.529411765,
    'rU:dC,6': 0.571428571,
    'rU:dG,14': 0.28571428600000004,
    'rU:dT,18': 0.666666667,
    'rA:dG,13': 0.21052631600000002,
    'rA:dG,12': 0.263157895,
    'rA:dG,11': 0.4,
    'rA:dG,10': 0.333333333,
    'rA:dA,19': 0.538461538,
    'rA:dA,18': 0.5,
    'rA:dG,15': 0.272727273,
    'rA:dG,14': 0.214285714,
    'rA:dA,15': 0.2,
    'rA:dA,14': 0.533333333,
    'rA:dA,17': 0.133333333,
    'rA:dA,16': 0.0,
    'rA:dA,11': 0.307692308,
    'rA:dA,10': 0.882352941,
    'rA:dA,13': 0.3,
    'rA:dA,12': 0.333333333,
    'rG:dA,13': 0.3,
    'rG:dA,12': 0.384615385,
    'rG:dA,11': 0.384615385,
    'rG:dA,10': 0.8125,
    'rG:dA,17': 0.25,
    'rG:dA,16': 0.0,
    'rG:dA,15': 0.14285714300000002,
    'rG:dA,6': 0.666666667,
    'rG:dG,20': 0.428571429,
    'rG:dA,19': 0.666666667,
    'rG:dA,18': 0.666666667,
    'rU:dC,4': 0.625,
    'rG:dT,12': 0.933333333,
    'rG:dT,13': 0.923076923,
    'rU:dG,11': 0.666666667,
    'rC:dA,3': 0.6875,
    'rC:dA,2': 0.9090909090000001,
    'rC:dA,1': 1.0,
    'rC:dA,7': 0.8125,
    'rC:dA,6': 0.9285714290000001,
    'rC:dA,5': 0.636363636,
    'rC:dA,4': 0.8,
    'rC:dA,9': 0.875,
    'rC:dA,8': 0.875,
    'rU:dT,6': 0.8666666670000001,
    'rA:dG,20': 0.22727272699999998,
    'rG:dT,18': 0.692307692,
    'rU:dG,10': 0.533333333,
    'rG:dT,19': 0.7142857140000001,
    'rG:dA,20': 0.7,
    'rC:dT,20': 0.5,
    'rU:dC,2': 0.84,
    'rG:dG,10': 0.4,
    'rC:dA,17': 0.46666666700000003,
    'rC:dA,16': 0.307692308,
    'rC:dA,15': 0.066666667,
    'rC:dA,14': 0.7333333329999999,
    'rC:dA,13': 0.7,
    'rC:dA,12': 0.538461538,
    'rC:dA,11': 0.307692308,
    'rC:dA,10': 0.9411764709999999,
    'rG:dG,11': 0.428571429,
    'rU:dC,20': 0.176470588,
    'rG:dG,3': 0.384615385,
    'rC:dA,19': 0.46153846200000004,
    'rC:dA,18': 0.642857143,
    'rU:dG,17': 0.705882353,
    'rU:dG,16': 0.666666667,
    'rU:dG,15': 0.272727273,
    'rG:dG,2': 0.692307692,
    'rU:dG,13': 0.7894736840000001,
    'rU:dG,12': 0.947368421,
    'rG:dA,9': 0.533333333,
    'rG:dA,8': 0.625,
    'rG:dA,7': 0.571428571,
    'rG:dG,5': 0.7857142859999999,
    'rG:dA,5': 0.3,
    'rG:dA,4': 0.363636364,
    'rG:dA,3': 0.5,
    'rG:dA,2': 0.636363636,
    'rG:dA,1': 1.0,
    'rG:dG,4': 0.529411765,
    'rG:dG,1': 0.7142857140000001,
    'rA:dC,9': 0.666666667,
    'rG:dG,7': 0.6875,
    'rG:dT,5': 0.8666666670000001,
    'rU:dT,20': 0.5625,
    'rC:dC,15': 0.05,
    'rC:dC,14': 0.0,
    'rC:dC,17': 0.058823529000000006,
    'rC:dC,16': 0.153846154,
    'rC:dC,11': 0.25,
    'rC:dC,10': 0.38888888899999996,
    'rC:dC,13': 0.13636363599999998,
    'rC:dC,12': 0.444444444,
    'rC:dA,20': 0.3,
    'rC:dC,19': 0.125,
    'rC:dC,18': 0.133333333,
    'rA:dA,1': 1.0,
    'rA:dA,3': 0.705882353,
    'rA:dA,2': 0.727272727,
    'rA:dA,5': 0.363636364,
    'rA:dA,4': 0.636363636,
    'rA:dA,7': 0.4375,
    'rA:dA,6': 0.7142857140000001,
    'rA:dA,9': 0.6,
    'rA:dA,8': 0.428571429,
    'rU:dG,20': 0.090909091,
    'rC:dC,9': 0.6190476189999999,
    'rC:dC,8': 0.642857143,
    'rU:dT,10': 0.857142857,
    'rU:dT,11': 0.75,
    'rU:dT,16': 0.9090909090000001,
    'rU:dT,17': 0.533333333,
    'rU:dT,14': 0.6190476189999999,
    'rU:dT,15': 0.578947368,
    'rC:dC,1': 0.913043478,
    'rU:dT,3': 0.7142857140000001,
    'rC:dC,3': 0.5,
    'rC:dC,2': 0.695652174,
    'rC:dC,5': 0.6,
    'rC:dC,4': 0.5,
    'rC:dC,7': 0.470588235,
    'rC:dC,6': 0.5,
    'rU:dT,4': 0.47619047600000003,
    'rU:dT,8': 0.8,
    'rU:dT,9': 0.9285714290000001,
    'rA:dC,19': 0.375,
    'rA:dC,18': 0.4,
    'rA:dC,17': 0.176470588,
    'rA:dC,16': 0.192307692,
    'rA:dC,15': 0.65,
    'rA:dC,14': 0.46666666700000003,
    'rA:dC,13': 0.6521739129999999,
    'rA:dC,12': 0.7222222220000001,
    'rA:dC,11': 0.65,
    'rA:dC,10': 0.5555555560000001,
    'rU:dC,7': 0.588235294,
    'rC:dT,8': 0.65,
    'rC:dT,9': 0.857142857,
    'rC:dT,6': 0.9285714290000001,
    'rC:dT,7': 0.75,
    'rC:dT,4': 0.842105263,
    'rC:dT,5': 0.571428571,
    'rC:dT,2': 0.727272727,
    'rC:dT,3': 0.8666666670000001,
    'rC:dT,1': 1.0,
    'rA:dC,8': 0.7333333329999999,
    'rU:dT,1': 1.0,
    'rU:dC,3': 0.5,
    'rU:dC,1': 0.956521739,
    'rU:dT,2': 0.846153846,
    'rU:dG,19': 0.275862069,
    'rG:dT,14': 0.75,
    'rG:dT,15': 0.9411764709999999,
    'rG:dT,16': 1.0,
    'rG:dT,17': 0.933333333,
    'rG:dT,10': 0.933333333,
    'rG:dT,11': 1.0,
    'rA:dG,9': 0.571428571,
    'rA:dG,8': 0.428571429,
    'rA:dG,7': 0.4375,
    'rA:dG,6': 0.454545455,
    'rA:dG,5': 0.5,
    'rA:dG,4': 0.352941176,
    'rA:dG,3': 0.428571429,
    'rA:dG,2': 0.7857142859999999,
    'rA:dG,1': 0.857142857,
    'rU:dT,5': 0.5,
    'rG:dT,2': 0.846153846,
    'rA:dC,3': 0.611111111,
    'rA:dC,20': 0.764705882,
    'rG:dT,1': 0.9,
    'rG:dT,6': 1.0,
    'rG:dT,7': 1.0,
    'rG:dT,4': 0.9,
    'rC:dT,19': 0.428571429,
    'rG:dG,9': 0.538461538,
    'rG:dG,8': 0.615384615,
    'rG:dT,8': 1.0,
    'rG:dT,9': 0.642857143,
    'rU:dG,18': 0.428571429,
    'rU:dT,7': 0.875,
    'rG:dG,6': 0.681818182,
    'rA:dA,20': 0.6,
    'rU:dC,9': 0.6190476189999999,
    'rA:dG,17': 0.176470588,
    'rU:dC,8': 0.7333333329999999,
    'rA:dG,16': 0.0,
    'rA:dG,19': 0.20689655199999998,
    'rG:dT,3': 0.75,
    'rU:dG,3': 0.428571429,
    'rU:dG,2': 0.857142857,
    'rU:dG,1': 0.857142857,
    'rA:dG,18': 0.19047619,
    'rU:dG,7': 0.6875,
    'rU:dG,6': 0.9090909090000001,
    'rU:dG,5': 1.0,
    'rU:dG,4': 0.647058824,
    'rU:dG,9': 0.923076923,
    'rU:dG,8': 1.0,
    'rU:dC,19': 0.25,
    'rU:dC,18': 0.333333333,
    'rU:dC,13': 0.260869565,
    'rU:dC,12': 0.5,
    'rU:dC,11': 0.4,
    'rU:dC,10': 0.5,
    'rU:dC,17': 0.117647059,
    'rU:dC,16': 0.346153846,
    'rU:dC,15': 0.05,
    'rU:dC,14': 0.0,
    'rC:dT,10': 0.8666666670000001,
    'rC:dT,11': 0.75,
    'rC:dT,12': 0.7142857140000001,
    'rC:dT,13': 0.384615385,
    'rC:dT,14': 0.35,
    'rC:dT,15': 0.222222222,
    'rC:dT,16': 1.0,
    'rC:dT,17': 0.46666666700000003,
    'rC:dT,18': 0.538461538,
    'rA:dC,2': 0.8,
    'rA:dC,1': 1.0,
    'rA:dC,7': 0.705882353,
    'rA:dC,6': 0.7142857140000001,
    'rA:dC,5': 0.72,
    'rA:dC,4': 0.625,
    'rU:dT,19': 0.28571428600000004,
}

PAM_SCORES = {
    'AA': 0.0,
    'AC': 0.0,
    'GT': 0.016129031999999998,
    'AG': 0.25925925899999996,
    'CC': 0.0,
    'CA': 0.0,
    'CG': 0.107142857,
    'TT': 0.0,
    'GG': 1.0,
    'GC': 0.022222222000000003,
    'AT': 0.0,
    'GA': 0.06944444400000001,
    'TG': 0.038961038999999996,
    'TA': 0.0,
    'TC': 0.0,
    'CT': 0.0,
}

if __name__ == "__main__":
    main()
