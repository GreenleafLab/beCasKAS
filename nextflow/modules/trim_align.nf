// Reproduces 1_trim_align.sh + 2_align_helper.sh from the legacy repo.
// Per-sample: trim_galore --paired --clip_R2 ${params.clip_r2} (with --fastqc for QC reports),
// then KAS-Analyzer KAS-seq. KAS-Analyzer dispatches to bowtie2 by default (no `-a` flag).

process TRIM_ALIGN {

    tag "${sample.id}"
    label 'process_high'

    input:
    tuple val(sample), path(r1), path(r2)
    path refdir                   // flat dir containing hg38.{1-4,rev.1,rev.2}.bt2 alongside hg38.fasta and other refs

    output:
    tuple val(sample), path("${sample.id}/Bam_files/${sample.id}_sorted.bam"),     emit: bam
    tuple val(sample), path("${sample.id}/Bam_files/${sample.id}_sorted.bam.bai"), emit: bai
    tuple val(sample), path("${sample.id}/Bed_files/${sample.id}.bed"),            emit: bed
    tuple val(sample), path("${sample.id}/BedGraph_files/${sample.id}.bg"),        emit: bg
    tuple val(sample), path("${sample.id}/Summary/*"),                           emit: summary, optional: true
    path("*_fastqc.{zip,html}"),                                             emit: fastqc, optional: true

    script:
    """
    set -euo pipefail

    # Decompress only if gzipped (KAS-Analyzer expects .fastq, not .fastq.gz)
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

    # Trim with trim_galore; --fastqc enables FastQC reports on both raw and trimmed reads
    # (matches original 2_align_helper.sh, which delegates fastqc to trim_galore via --fastqc, not separate KAS-Analyzer fastqc calls)
    trim_galore -j ${task.cpus} --fastqc --paired --clip_R2 ${params.clip_r2} \$R1 \$R2

    R1_TRIM=\$(basename \$R1 .fastq)_val_1.fq
    R2_TRIM=\$(basename \$R2 .fastq)_val_2.fq

    # KAS-seq align: -i takes the bowtie2 INDEX PREFIX (not a directory); -s is the assembly id literal "hg38"
    REFDIR_ABS=\$(readlink -f ${refdir})
    KAS-Analyzer KAS-seq -t ${task.cpus} -i \$REFDIR_ABS/hg38 -o ${sample.id} -s hg38 -1 \$R1_TRIM -2 \$R2_TRIM
    """
}
