// The broad arm's interval set, built with the same consensus rule DiffBind applies to the sharp
// peaks: merge the peaksets of every sample in this case's contrast (treatment replicates plus the
// shared controls) and keep intervals present in >= min_overlap of them. min_overlap is the
// effective per-case minOverlap DiffBind's dba.count used for this case, supplied by the caller,
// so with real replicates the differential peaks and these broad regions share a threshold and a
// peakset membership. That is the consensus rule only, not the whole of dba.count: its default
// read-count filter is not mirrored here. Membership also diverges under pseudo-replicates:
// when a condition has n=1,
// build_diffbind_samplesheet.py gives every pseudo-rep the parent's peakset, so DiffBind counts
// that one peakset once per pseudo-rep while this module counts it once.
//
// Deliberately NOT matched: dba.count's summits = 200, which re-centres every consensus peak
// to a 401-bp window. That is right for sharp peaks and would destroy a broad domain. The sharp
// non-differential arm does not summit-recentre either.
//
// Nothing is subtracted: the broad arm runs over the whole consensus, independent of the sharp
// arms, so a region can be scored by more than one arm. ARM_OVERLAP reports how much they share.

process BROAD_REGIONS {

    tag "${case_id}"

    input:
    tuple val(case_id), path(broad_peaks, stageAs: 'broad_??.bed'), val(min_overlap)

    output:
    tuple val(case_id), path("broad_regions.bed"), emit: regions
    tuple val(case_id), path("broad_regions_summary.txt"), emit: summary

    script:
    def peaks      = [broad_peaks].flatten()
    def n_peaksets = peaks.size()
    """
    set -euo pipefail

    if [ ${min_overlap} -gt ${n_peaksets} ]; then
        echo "ERROR: --min_overlap=${min_overlap} exceeds the ${n_peaksets} peakset(s) for case ${case_id}; no interval can satisfy it." >&2
        exit 1
    fi

    # EPIC2 BEDs carry a '#Chromosome...' header; strip it and tag each interval with its peakset
    i=0
    for f in ${peaks.join(' ')}; do
        i=\$((i + 1))
        awk -v r=\$i 'BEGIN{OFS="\\t"} !/^#/ && NF >= 3 { print \$1, \$2, \$3, r }' "\$f"
    done | sort -k1,1 -k2,2n > tagged.bed

    # Consensus, same rule as dba.count(minOverlap): present in >= min_overlap of the peaksets.
    # merge is transitive: a chain of overlaps counts every peakset in it, even where no single
    # base lies in that many (as in DiffBind's own merge)
    if [ -s tagged.bed ]; then
        bedtools merge -i tagged.bed -c 4 -o count_distinct \\
            | awk -v m=${min_overlap} 'BEGIN{OFS="\\t"} \$4 >= m { print \$1, \$2, \$3 }' \\
            > broad_regions.bed
    else
        : > broad_regions.bed
    fi

    bp() { awk '{ bp += \$3 - \$2; n++ } END { printf "%d intervals, %.2f Mb", n, bp/1e6 }' "\$1"; }
    {
        echo "# broad arm regions: ${case_id}"
        echo "peaksets (trt + control):        ${n_peaksets}"
        echo "minOverlap (as DiffBind):        ${min_overlap}"
        echo "broad consensus (this file):     \$(bp broad_regions.bed)"
    } > broad_regions_summary.txt
    cat broad_regions_summary.txt
    """
}

// Each arm's variant-calling intervals (differential peaks; each non-differential arm's scored
// peaks), tabulated and drawn as a Venn over every combination of arms. sharp_bed is [] for a
// case the sharp arm skipped (an empty circle); without --sharp_nondiff there is no sharp circle.
process ARM_OVERLAP {

    tag "${case_id}"

    input:
    tuple val(case_id), path(diff_bed, stageAs: 'differential.bed'), path(sharp_bed, stageAs: 'sharp_nondiff.bed'), path(broad_bed, stageAs: 'broad_nondiff.bed')

    output:
    tuple val(case_id), path("arm_overlap.tsv"),      emit: table
    tuple val(case_id), path("arm_overlap_venn.pdf"), emit: venn

    script:
    def sharp = params.sharp_nondiff ? "sharp_nondiff=${sharp_bed ?: '/dev/null'}" : ''
    """
    ${projectDir}/nextflow/bin/arm_overlap.R ${case_id} differential=${diff_bed} ${sharp} broad_nondiff=${broad_bed}
    """
}
