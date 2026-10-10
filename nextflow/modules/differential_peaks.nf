// Cohort-level differential peak analysis (DiffBind, DESeq2 backend).

process BUILD_DIFFBIND_SAMPLESHEET {
    input:
    val rows         // list of [sample_id, condition, replicate, bam_path]
    val peaks_map    // list of [sample_id, peaks_path] 
    output:
    path 'diffbind_samplesheet_*.csv', emit: samplesheets

    script:
    """
    set -euo pipefail

    cat > rows.json <<'JSON'
    ${groovy.json.JsonOutput.toJson(rows)}
    JSON
    cat > peaks.json <<'JSON'
    ${groovy.json.JsonOutput.toJson(peaks_map)}
    JSON

    ${projectDir}/nextflow/bin/build_diffbind_samplesheet.py \\
        --rows rows.json \\
        --peaks peaks.json \\
        --output-dir .
    """
}

process DIFFERENTIAL_PEAKS {

    label 'process_high'
    tag "${case_id}"

    input:
    tuple val(case_id), path(samplesheet), val(params_json)

    output:
    tuple val(case_id), path("differential_peaks.bed"),      emit: bed
    tuple val(case_id), path("nondiff_peaks.bed"),           emit: nondiff_peaks
    tuple val(case_id), path("differential_peaks_full.tsv"), emit: full
    tuple val(case_id), path("differential_peaks_method.txt"), emit: method
    tuple val(case_id), path("volcano.pdf"),                 emit: volcano
    tuple val(case_id), path("ma_plot.pdf"),                 emit: ma_plot
    tuple val(case_id), path("dba_object.rds"),              emit: dba_object

    script:
    """
    set -euo pipefail

    cat > diffbind_params.json <<'JSON'
    ${params_json}
    JSON

    ${projectDir}/nextflow/bin/differential_peaks.R \\
        --samplesheet ${samplesheet} \\
        --output-dir . \\
        --min-overlap ${params.min_overlap} \\
        --seed ${params.seed} \\
        --threads ${task.cpus} \\
        --n-pseudo-replicates ${params.n_pseudo_replicates} \\
        --diffbind-params diffbind_params.json
    """
}
