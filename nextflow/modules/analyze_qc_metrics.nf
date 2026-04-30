// Cohort QC: FRiP via KAS-Analyzer + correlation via deeptools multiBigwigSummary.
// Operates on non-input samples only (matches 3_analyze.sh:s1, s2 array).

process QC_METRICS {

    input:
    path(beds)         // non-input sample BEDs
    path(peaks)        // non-input sample peak BEDs (from PEAK_CALL)
    path(bigwigs)      // non-input sample .nor.bigWig (from NORMALIZE_BIGWIG)
    val sample_labels  // ordered list matching beds/peaks/bigwigs

    output:
    path("KAS-Seq_FRiP.txt"),       emit: frip,        optional: true
    path("all-deeptools.pdf"),      emit: correlation_pdf
    path("all.tab"),                emit: correlation_tab
    path("all.npz"),                emit: correlation_npz

    script:
    def label_str = sample_labels.join('\\n')
    """
    set -euo pipefail

    # Build FRiP input lists
    printf '%s\\n' ${beds.collect{"$it"}.join(' ')} | sed 's| |\\n|g' > KAS-seq_input.txt
    printf '${label_str}' > labels_input.txt
    printf '%s\\n' ${peaks.collect{"$it"}.join(' ')} | sed 's| |\\n|g' > peaks_input_files.txt

    KAS-Analyzer FRiP -o KAS-Seq -p peaks_input_files.txt -l labels_input.txt -k KAS-seq_input.txt || true

    # Correlation matrix (deeptools)
    multiBigwigSummary bins -b ${bigwigs.join(' ')} -o all.npz --outRawCounts all.tab

    plotCorrelation -in all.npz -c pearson -p scatterplot -o all-deeptools.pdf \\
        --removeOutliers --log1p --xRange 0 4 --yRange 0 4
    """
}
