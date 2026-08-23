// Most-likely gRNA binding site per differential peak, via GuideScan2 + Doench 2016 CFD.

process GUIDESCAN_BINDING_SITES {

    tag "${case_id}"

    input:
    tuple val(case_id), path(differential_peaks_bed)
    path fasta
    path fasta_fai

    output:
    tuple val(case_id), path("binding_sites.tsv"), emit: binding_sites

    script:
    """
    set -euo pipefail

    ${projectDir}/nextflow/bin/guidescan.py \\
        --bed ${differential_peaks_bed} \\
        --genome ${fasta} \\
        --grna ${params.guidescan_grna} --pam ${params.guidescan_pam} \\
        --mismatches ${params.guidescan_mismatches} \\
        --rna-bulges ${params.guidescan_rna_bulges} --dna-bulges ${params.guidescan_dna_bulges} \\
        --chunk-size ${params.guidescan_chunk_size} --workers ${task.cpus} \\
        --outfile binding_sites.tsv
    """
}