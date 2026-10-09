// Most-likely gRNA binding site per peak (Doench 2016 CFD), in two stages.
// Stage 1 anchors on an exact seed match and scores that one placement
// Stage 2 runs GuideScan2 enumerate over remaining peaks with no exact seed anchor
//
// Both the differential arm and the non-differential arm (see main.nf) run these processes; the
// arms differ only in the peak BED they are given and in whether stage 2 runs at all.

// Per-case gRNA/seed settings: cases in one run can target different loci, so each
// field comes from --guidescan_params[case_id] if set, else the global --guidescan_*.
// Returns null when no gRNA is configured for this case.
def resolveGuidescanParams(Map overrides, String case_id) {
    def o    = (overrides[case_id] ?: [:])
    def grna = (o.grna ?: params.guidescan_grna)
    if (!grna) return null
    def pam  = (o.pam ?: params.guidescan_pam)

    // Default seed is a substring of spacer+PAM (PAM-proximal 5nt of spacer + PAM)
    def seed = (o.seed ?: params.guidescan_seed)
    if (!seed) {
        if (grna.toString().length() < 5) {
            error("guidescan: case '${case_id}' spacer '${grna}' is shorter than 5 nt, so no " +
                  "seed can be derived from it -- set one explicitly via --guidescan_seed.")
        }
        seed = grna.toString()[-5..-1] + pam
    }

    [ grna: grna, pam: pam, seed: seed,
      mismatches          : (o.mismatches           != null ? o.mismatches           : params.guidescan_mismatches),
      rna_bulges          : (o.rna_bulges           != null ? o.rna_bulges           : params.guidescan_rna_bulges),
      dna_bulges          : (o.dna_bulges           != null ? o.dna_bulges           : params.guidescan_dna_bulges),
      chunk_size          : (o.chunk_size           != null ? o.chunk_size           : params.guidescan_chunk_size),
      summit_width        : (o.summit_width         != null ? o.summit_width         : params.guidescan_summit_width),
      enumerate_mismatches: (o.enumerate_mismatches != null ? o.enumerate_mismatches : params.guidescan_enumerate_mismatches),
      enumerate_chunk_size: (o.enumerate_chunk_size != null ? o.enumerate_chunk_size : params.guidescan_enumerate_chunk_size) ]
}

process GUIDESCAN_SEED {

    tag "${case_id}"

    input:
    tuple val(case_id), path(peaks_bed), val(gs)
    path fasta
    path fasta_fai

    output:
    tuple val(case_id), path("binding_sites_seed.tsv"),          emit: binding_sites
    tuple val(case_id), path("binding_sites_seed_raw_hits.tsv"), emit: raw_hits
    tuple val(case_id), path("no_seed_peaks.bed"), val(gs),      emit: no_seed_bed
    tuple val(case_id), path("scored_peaks.bed"),                emit: scored_peaks

    script:
    """
    set -euo pipefail

    ${projectDir}/nextflow/bin/guidescan_seed.py \\
        --bed ${peaks_bed} \\
        --genome ${fasta} \\
        --grna ${gs.grna} --pam ${gs.pam} --seed ${gs.seed} \\
        --mismatches ${gs.mismatches} \\
        --rna-bulges ${gs.rna_bulges} --dna-bulges ${gs.dna_bulges} \\
        --chunk-size ${gs.chunk_size} --workers ${task.cpus} \\
        --summit-width ${gs.summit_width} \\
        --workdir . \\
        --outfile binding_sites_seed.tsv \\
        --raw-hits-outfile binding_sites_seed_raw_hits.tsv

    # Stage 2's input: peaks with no exact seed anchor anywhere, on either strand.
    # seed_hits is col 11; the peak's own coordinates are cols 12-14.
    awk -F'\\t' 'NR > 1 && \$11 == 0 { print \$12 "\\t" \$13 "\\t" \$14 }' \\
        binding_sites_seed.tsv | sort -k1,1 -k2,2n -u > no_seed_peaks.bed

    # Build interval list for non-differential arm (drop unscored sites)
    awk -F'\\t' 'NR > 1 && \$2 != "." { print \$12 "\\t" \$13 "\\t" \$14 }' \\
        binding_sites_seed.tsv | sort -k1,1 -k2,2n -u > scored_peaks.bed
    """
}

process GUIDESCAN_ENUMERATE {

    tag "${case_id}"

    input:
    tuple val(case_id), path(no_seed_bed), val(gs)
    path fasta
    path fasta_fai

    output:
    tuple val(case_id), path("binding_sites_enumerate.tsv"), emit: binding_sites

    script:
    """
    set -euo pipefail

    ${projectDir}/nextflow/bin/guidescan.py \\
        --bed ${no_seed_bed} \\
        --genome ${fasta} \\
        --grna ${gs.grna} --pam ${gs.pam} \\
        --mismatches ${gs.enumerate_mismatches} \\
        --rna-bulges ${gs.rna_bulges} --dna-bulges ${gs.dna_bulges} \\
        --chunk-size ${gs.enumerate_chunk_size} \\
        --workers ${task.cpus} --threads 1 \\
        --workdir . \\
        --outfile binding_sites_enumerate.tsv
    """
}