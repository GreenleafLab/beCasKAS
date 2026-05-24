// Combined module for 3_analyze.sh: per-sample peak calling, cohort RPM statistics,
// normalize + bigWig tracks, and cohort QC (FRiP + correlation). Each process is
// independently importable from main.nf:
//   include { PEAK_CALL; KAS_STATISTICS; NORMALIZE_BIGWIG; QC_METRICS } from './nextflow/modules/analyze'

// =====================================================================================
// PEAK_CALL -- per-sample MACS2 peak calling via KAS-Analyzer
// Reproduces 3_analyze.sh:55-58 (peakscalling + post-process to BED).
// =====================================================================================

process PEAK_CALL {

    tag "${sample.id}"
    label 'process_medium'

    input:
    tuple val(sample), path(sample_bed), path(control_bed)
    val genome_name

    output:
    tuple val(sample), path("${sample.id}_input_macs2_peaks.bed"),          emit: peaks
    tuple val(sample), path("${sample.id}_input_macs2_peaks_abridged.bed"), emit: peaks_abridged
    tuple val(sample), path("${sample.id}_input_peaks.xls"),                emit: macs2_xls

    script:
    """
    set -euo pipefail

    KAS-Analyzer peakscalling \\
        -t ${sample_bed} \\
        -c ${control_bed} \\
        -b ${params.macs2_mode} \\
        -o ${sample.id}_input \\
        -g ${genome_name}

    # Strip MACS2 xls comment lines + blank rows -> headered BED (used by DiffBind)
    awk '!/^#/ && NF' ${sample.id}_input_peaks.xls > ${sample.id}_input_macs2_peaks.bed

    # Strip the column-header line -> pure BED (used by FRiP / bedtools)
    awk 'NR > 1' ${sample.id}_input_macs2_peaks.bed > ${sample.id}_input_macs2_peaks_abridged.bed
    """
}


// =====================================================================================
// KAS_STATISTICS -- cohort RPM ratio computation
// Reproduces 3_analyze.sh:74-96 verbatim:
//   1. KAS-Analyzer statistics aggregates per-sample summaries into a cohort summary
//   2. awk #1: appends RPM column = 1,000,000 / Deduplicated_reads (column 4)
//   3. awk #2: extracts the RPM column for samples -> ratios.txt
//
// The only generalization vs the original: the row-count boundary in the second awk
// changes from the hardcoded "NR <= 16" (which assumed 1 header + 15 samples in the
// preprint's specific cohort) to "NR <= n_samples + 1" so cohorts of any size work.
// The RPM math is bit-identical.
// =====================================================================================

process KAS_STATISTICS {

    label 'process_low'

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
    mkdir -p Summary
    for f in ${summary_files}; do
        cp "\$f" Summary/
    done

    # Cohort-level statistics -> produces Summary/KAS_statistics_summary.txt (verbatim from 3_analyze.sh:74)
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


// =====================================================================================
// NORMALIZE_BIGWIG -- RPM-ratio normalization + UCSC bigWig tracks
// Reproduces the normalize + UCSC subblock of 3_analyze.sh (lines 132-138).
// Cohort-aware per-sample: needs all samples' bg files + a ratios.txt from KAS_STATISTICS.
// =====================================================================================

process NORMALIZE_BIGWIG {

    label 'process_medium'

    input:
    path(bg_files)         // collected: all samples' .bg
    path(ratios_file)      // ratios.txt with one RPM per sample line, in the same order as bg_files
    path(labels_file)      // UCSC track labels
    val genome_name

    output:
    path("*.nor.bigWig"),         emit: bigwigs
    path("*.nor.bg"),             emit: nor_bg
    path("KAS-seq.txt"),          emit: bg_list
    path("UCSC_track_names.txt"), emit: ucsc_labels

    script:
    """
    set -euo pipefail

    # Build the file lists KAS-Analyzer expects (one filename per line).
    ls *.bg | grep -v '\\.nor\\.bg\$' > KAS-seq.txt
    [ "${labels_file}" != "UCSC_track_names.txt" ] && cp ${labels_file} UCSC_track_names.txt || true

    # Normalize via KAS-Analyzer with provided ratios.txt
    KAS-Analyzer normalize -k KAS-seq.txt -r ${ratios_file} -s ${genome_name} -b

    # Track names file for UCSC
    ls *.nor.bg > KAS-seq-nor.txt
    KAS-Analyzer UCSC -k KAS-seq-nor.txt -l UCSC_track_names.txt
    """
}


// =====================================================================================
// QC_METRICS -- cohort FRiP + multiBigwigSummary correlation
// Reproduces 3_analyze.sh:102 (FRiP) and 144-150 (deeptools correlation).
// Operates on non-input samples only (matches the s1, s2 array filter in the original).
// =====================================================================================

process QC_METRICS {

    label 'process_medium'

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
