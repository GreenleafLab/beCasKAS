// Per-sample dedup: stamp read groups, mark duplicates with picard, then collapse to consensus reads with gencore.

process DEDUP_CONSENSUS {

    tag "${sample.id}"
    label 'process_medium'

    input:
    tuple val(sample), path(sorted_bam)
    path fasta
    path fasta_fai

    output:
    tuple val(sample), path("${sample.id}_picard-gencore-dedup.bam"),     emit: bam
    tuple val(sample), path("${sample.id}_picard-gencore-dedup.bam.bai"), emit: bai
    tuple val(sample), path("${sample.id}_picard-gencore-dedup.html"),    emit: html
    tuple val(sample), path("${sample.id}_marked_metrics.txt"),           emit: dup_metrics,  optional: true
    tuple val(sample), path("gencore.json"),                              emit: gencore_json, optional: true

    script:
    """
    set -euo pipefail

    # Add default read group (required for downstream picard)
    picard AddOrReplaceReadGroups \\
        I=${sorted_bam} \\
        O=${sample.id}_sorted_rg.bam \\
        RGID=4 \\
        RGLB=lib1 \\
        RGPL=ILLUMINA \\
        RGPU=unit1 \\
        RGSM=20

    # Mark duplicates
    picard MarkDuplicates \\
        I=${sample.id}_sorted_rg.bam \\
        O=${sample.id}_marked.bam \\
        M=${sample.id}_marked_metrics.txt

    # Consensus dedup via gencore
    gencore -i ${sample.id}_marked.bam \\
        -o ${sample.id}_picard-gencore-dedup.bam \\
        -r ${fasta} \\
        --html ${sample.id}_picard-gencore-dedup.html

    samtools index -@ ${task.cpus} ${sample.id}_picard-gencore-dedup.bam ${sample.id}_picard-gencore-dedup.bam.bai

    rm ${sample.id}_sorted_rg.bam ${sample.id}_marked.bam
    """
}
