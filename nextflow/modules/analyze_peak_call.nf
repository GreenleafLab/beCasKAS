// Reproduces the peak-calling block of 3_analyze.sh.

process PEAK_CALL {

    tag "${meta.id}"

    input:
    tuple val(meta), path(sample_bed), path(control_bed)
    val genome_name

    output:
    tuple val(meta), path("${meta.id}_input_macs2_peaks.bed"),          emit: peaks
    tuple val(meta), path("${meta.id}_input_macs2_peaks_abridged.bed"), emit: peaks_abridged
    tuple val(meta), path("${meta.id}_input_peaks.xls"),                emit: macs2_xls

    script:
    """
    set -euo pipefail

    KAS-Analyzer peakscalling \\
        -t ${sample_bed} \\
        -c ${control_bed} \\
        -b ${params.macs2_mode} \\
        -o ${meta.id}_input \\
        -g ${genome_name}

    # Strip MACS2 xls comment lines + blank rows → headered BED (used by DiffBind)
    awk '!/^#/ && NF' ${meta.id}_input_peaks.xls > ${meta.id}_input_macs2_peaks.bed

    # Strip the column-header line → pure BED (used by FRiP / bedtools)
    awk 'NR > 1' ${meta.id}_input_macs2_peaks.bed > ${meta.id}_input_macs2_peaks_abridged.bed
    """
}
