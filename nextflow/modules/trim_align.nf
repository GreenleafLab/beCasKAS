// Per-sample read trimming + bowtie2 alignment:
// trim_galore --paired --clip_R2 with FastQC on raw + trimmed reads, then KAS-Analyzer KAS-seq (bowtie2 backend).

process TRIM_ALIGN {

    tag "${sample.id}"
    label 'process_high'

    input:
    tuple val(sample), path(r1), path(r2)
    path refdir                   // flat dir containing hg38.{1-4,rev.1,rev.2}.bt2 alongside hg38.fasta and other refs

    output:
    tuple val(sample), path("Bam_files/${sample.id}_sorted.bam"),     emit: bam
    tuple val(sample), path("Bam_files/${sample.id}_sorted.bam.bai"), emit: bai
    tuple val(sample), path("Bed_files/${sample.id}.bed"),            emit: bed
    tuple val(sample), path("BedGraph_files/${sample.id}.bg"),        emit: bg
    tuple val(sample), path("Summary/*"),                             emit: summary, optional: true
    path("*_fastqc.{zip,html}"),                                      emit: fastqc, optional: true

    script:
    """
    set -euo pipefail

    # Decompress (for KAS-Analyzer))
    R1=${r1}
    R2=${r2}
    if [[ "\$R1" == *.gz ]]; then
        gunzip -c \$R1 > ${sample.id}_R1.fastq
        R1=${sample.id}_R1.fastq
    fi
    if [[ "\$R2" == *.gz ]]; then
        gunzip -c \$R2 > ${sample.id}_R2.fastq
        R2=${sample.id}_R2.fastq
    fi

    # Trim with trim_galore
    trim_galore -j ${task.cpus} --fastqc --paired --clip_R2 ${params.clip_r2} \$R1 \$R2

    R1_TRIM=\$(basename \$R1 .fastq)_val_1.fq
    R2_TRIM=\$(basename \$R2 .fastq)_val_2.fq

    # KAS-seq align"
    REFDIR_ABS=\$(readlink -f ${refdir})
    KAS-Analyzer KAS-seq -t ${task.cpus} -i \$REFDIR_ABS/hg38 -o ${sample.id} -s hg38 -1 \$R1_TRIM -2 \$R2_TRIM
    """
}
