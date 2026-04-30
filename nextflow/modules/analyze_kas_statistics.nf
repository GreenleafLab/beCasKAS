// Reproduces 3_analyze.sh:74-96 verbatim:
//   1. KAS-Analyzer statistics aggregates per-sample summaries into a cohort summary
//   2. awk #1: appends RPM column = 1,000,000 / Deduplicated_reads (column 4)
//   3. awk #2: extracts the RPM column for samples → ratios.txt
//
// The only generalization vs the original: the row-count boundary in the second awk
// changes from the hardcoded "NR <= 16" (which assumed 1 header + 15 samples in the
// preprint's specific cohort) to "NR <= n_samples + 1" so cohorts of any size work.
// The RPM math is bit-identical.

process KAS_STATISTICS {

    input:
    path summary_files   // flat list of all per-sample summary files emitted by TRIM_ALIGN
    val sample_labels    // ordered list of sample IDs (one per sample, matches FRiP/UCSC label order)

    output:
    path 'ratios.txt',                          emit: ratios
    path 'KAS_statistics_summary.txt',          emit: stats
    path 'KAS_statistics_summary_rpm.txt',      emit: stats_rpm
    path 'UCSC_track_names.txt',                emit: labels

    script:
    def n_samples = sample_labels.size()
    """
    set -euo pipefail

    # Stage all per-sample summary files into a single Summary/ dir as KAS-Analyzer expects.
    # Assumes KAS-Analyzer names summary files with sample-specific prefixes (it does, e.g.
    # ${sample_id}_*); if collisions occur the impl should rename per-sample at staging.
    mkdir -p Summary
    for f in ${summary_files}; do
        cp "\$f" Summary/
    done

    # Cohort-level statistics → produces Summary/KAS_statistics_summary.txt (verbatim from 3_analyze.sh:74)
    KAS-Analyzer statistics -o KAS_statistics -s ./Summary

    # awk #1 verbatim from 3_analyze.sh:79-91
    awk 'BEGIN { FS=OFS="\\t" }
        NR==1 {
            print \$0, "RPM"   # Print the header with the new column "RPM"
        }
        NR > 1 {
            if (\$4 > 0) {  # needs to know correct column
                rpm = 1000000 / \$4
            } else {
                rpm = 0  # If Deduplicated_reads is zero or negative, set RPM to 0
            }
            print \$0, rpm  # Print the entire row with the new RPM value
        }
    ' Summary/KAS_statistics_summary.txt > Summary/KAS_statistics_summary_rpm.txt

    # awk #2 from 3_analyze.sh:96, generalized to cohort size
    awk -F"\\t" 'NR > 1 && NR <= ${n_samples + 1} {print \$7}' \\
        Summary/KAS_statistics_summary_rpm.txt > ratios.txt

    # Surface the summary files at process root for output capture
    cp Summary/KAS_statistics_summary.txt     KAS_statistics_summary.txt
    cp Summary/KAS_statistics_summary_rpm.txt KAS_statistics_summary_rpm.txt

    # Labels file (sample IDs in matching order to ratios.txt rows)
    printf '%s\\n' ${sample_labels.collect{ "\"${it}\"" }.join(' ')} > UCSC_track_names.txt
    """
}
