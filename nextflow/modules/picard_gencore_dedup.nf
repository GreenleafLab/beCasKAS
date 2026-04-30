// Reproduces 4_picard-gencore.sh per-sample.

process DEDUP_CONSENSUS {

    tag "${meta.id}"

    input:
    tuple val(meta), path(sorted_bam)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.id}_picard-gencore-dedup.bam"),     emit: bam
    tuple val(meta), path("${meta.id}_picard-gencore-dedup.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}_picard-gencore-dedup.html"),    emit: html

    script:
    """
    set -euo pipefail

    # Add default read group (picard requires it for downstream tools)
    picard AddOrReplaceReadGroups \\
        I=${sorted_bam} \\
        O=${meta.id}_sorted_rg.bam \\
        RGID=4 \\
        RGLB=lib1 \\
        RGPL=ILLUMINA \\
        RGPU=unit1 \\
        RGSM=20

    # Mark duplicates
    picard MarkDuplicates \\
        I=${meta.id}_sorted_rg.bam \\
        O=${meta.id}_marked.bam \\
        M=${meta.id}_marked_metrics.txt

    # Consensus dedup via gencore
    gencore -i ${meta.id}_marked.bam \\
        -o ${meta.id}_picard-gencore-dedup.bam \\
        -r ${fasta} \\
        --html ${meta.id}_picard-gencore-dedup.html

    samtools index ${meta.id}_picard-gencore-dedup.bam ${meta.id}_picard-gencore-dedup.bam.bai

    rm ${meta.id}_sorted_rg.bam ${meta.id}_marked.bam
    """
}
