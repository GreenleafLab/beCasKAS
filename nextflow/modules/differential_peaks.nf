// Cohort-level differential peak analysis (DiffBind, DESeq2 backend).

process BUILD_DIFFBIND_SAMPLESHEET {

    input:
    val rows         // list of [sample_id, condition, replicate, bam_path]
    val peaks_map    // list of [sample_id, peaks_path] (real samples only; pseudo-reps reuse parent peaks)

    output:
    path 'diffbind_samplesheet.csv'

    script:
    """
    python3 - <<'PYEOF' > diffbind_samplesheet.csv
    import csv, sys

    rows = ${groovy.json.JsonOutput.toJson(rows)}
    peaks_map = dict([(p[0], p[1]) for p in ${groovy.json.JsonOutput.toJson(peaks_map)}])

    writer = csv.writer(sys.stdout)
    writer.writerow(['sample_id', 'condition', 'replicate', 'bam_path', 'peaks_path'])
    for r in rows:
        sid, cond, rep, bam = r
        parent = sid.split('_pseudo')[0]
        peaks = peaks_map.get(sid) or peaks_map.get(parent) or ''
        writer.writerow([sid, cond, rep, bam, peaks])
    PYEOF
    """
}

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
