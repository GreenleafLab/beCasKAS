// Only invoked when a condition has n=1 sample. Partitions BAM into params.n_pseudo_replicates buckets.

process SPLIT_BAM_PSEUDOREPS {

    tag "${sample.id}"

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample.id}_pseudo*.bam"), path("${sample.id}_pseudo*.bam.bai"), emit: pseudoreps

    script:
    """
    set -euo pipefail

    ${projectDir}/nextflow/bin/split_bam_pseudoreps.py \\
        --input ${bam} \\
        --output-prefix ${sample.id} \\
        --n-buckets ${params.n_pseudo_replicates} \\
        --seed ${params.seed}
    """
}
