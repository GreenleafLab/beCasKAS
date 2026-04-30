// Only invoked when a condition has n=1 sample. Partitions BAM into params.n_pseudo_replicates buckets.

process SPLIT_BAM_PSEUDOREPS {

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_pseudo*.bam"), path("${meta.id}_pseudo*.bam.bai"), emit: pseudoreps

    script:
    """
    set -euo pipefail

    ${projectDir}/nextflow/bin/split_bam_pseudoreps.py \\
        --input ${bam} \\
        --output-prefix ${meta.id} \\
        --n-buckets ${params.n_pseudo_replicates} \\
        --seed ${params.seed}
    """
}
