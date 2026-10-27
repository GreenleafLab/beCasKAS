// Per-sample read trimming + bowtie2 alignment:
// trim_galore --paired --clip_R2 with FastQC on raw + trimmed reads, then KAS-Analyzer KAS-seq (bowtie2 backend).

process TRIM_ALIGN {

    tag "${sample.id}"
    label 'process_high'

    input:
    tuple val(sample), path(r1), path(r2)
    path refdir
    val assembly

    output:
    tuple val(sample), path("Bam_files/${sample.id}_sorted.bam"),     emit: bam
    tuple val(sample), path("Bam_files/${sample.id}_sorted.bam.bai"), emit: bai
    tuple val(sample), path("Bed_files/${sample.id}.bed"),            emit: bed
    tuple val(sample), path("BedGraph_files/${sample.id}.bg"),        emit: bg
    tuple val(sample), path("Summary/*"),                             emit: summary, optional: true
    path("*_fastqc.{zip,html}"),                                      emit: fastqc, optional: true
    path("*_trimming_report.txt"),                                    emit: trim_report, optional: true

    script:
    """
    set -euo pipefail

    # Trim with trim_galore. Both trim_galore (0.6.10) and KAS-Analyzer/KAS-pipe2
    trim_galore -j ${task.cpus} --fastqc --paired --clip_R2 ${params.clip_r2} --basename ${sample.id} ${r1} ${r2}
    R1_TRIM=\$(ls ${sample.id}_val_1.fq* | head -n1)
    R2_TRIM=\$(ls ${sample.id}_val_2.fq* | head -n1)

    # KAS-seq align
    REFDIR_ABS=\$(readlink -f ${refdir})
    KAS-Analyzer KAS-seq -t ${task.cpus} -i \$REFDIR_ABS/${params.genome} -o ${sample.id} -s ${assembly} -1 \$R1_TRIM -2 \$R2_TRIM
    rm -f ${sample.id}_val_*.fq* ${sample.id}/${sample.id}_val_*.fq* 2>/dev/null || true
    """
}
