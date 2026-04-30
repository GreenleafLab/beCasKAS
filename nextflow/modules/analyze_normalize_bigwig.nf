// Reproduces the normalize + UCSC subblock of 3_analyze.sh.
// Cohort-aware per-sample: needs all samples' bg files + a ratios.txt computed from cohort stats.

process NORMALIZE_BIGWIG {

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
    cp ${labels_file} UCSC_track_names.txt

    # Normalize via KAS-Analyzer with provided ratios.txt
    KAS-Analyzer normalize -m ratios -k KAS-seq.txt -r ${ratios_file} -s ${genome_name} -b

    # Track names file for UCSC
    ls *.nor.bg > KAS-seq-nor.txt
    KAS-Analyzer UCSC -k KAS-seq-nor.txt -l UCSC_track_names.txt
    """
}
