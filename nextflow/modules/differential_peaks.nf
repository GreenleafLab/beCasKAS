// Cohort-level differential peak analysis (DiffBind, DESeq2 backend).

process DIFFERENTIAL_PEAKS {

    label 'process_high'

    input:
    path samplesheet      // pipeline samplesheet with absolute bam_path, peaks_path columns added

    output:
    path("differential_peaks.bed"),         emit: bed
    path("differential_peaks_full.tsv"),    emit: full
    path("differential_peaks_method.txt"),  emit: method
    path("volcano.pdf"),                    emit: volcano
    path("dba_object.rds"),                 emit: dba_object

    script:
    """
    set -euo pipefail

    ${projectDir}/nextflow/bin/differential_peaks.R \\
        --samplesheet ${samplesheet} \\
        --output-dir . \\
        --fdr ${params.fdr} \\
        --lfc ${params.lfc} \\
        --min-overlap ${params.min_overlap} \\
        ${params.greylist ? '--greylist' : ''} \\
        --seed ${params.seed} \\
        --threads ${task.cpus} \\
        --n-pseudo-replicates ${params.n_pseudo_replicates}
    """
}
