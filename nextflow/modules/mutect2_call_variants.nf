// Reproduces 6_gatk-mutect2.sh per non-input sample.

process CALL_VARIANTS {

    tag "${sample.id}"
    label 'process_high'

    input:
    tuple val(sample), path(dedup_bam), path(dedup_bai)
    path differential_peaks_bed
    path fasta
    path fasta_fai
    path fasta_dict
    path gnomad
    path gnomad_tbi
    path pon
    path pon_tbi

    output:
    tuple val(sample), path("${sample.id}_no-alt-chr.vcf"), emit: vcf_no_alt
    tuple val(sample), path("${sample.id}_CT.vcf"),         emit: vcf_ct
    tuple val(sample), path("${sample.id}_AG.vcf"),         emit: vcf_ag
    tuple val(sample), path("${sample.id}.vcf"),            emit: vcf_raw

    script:
    """
    set -euo pipefail

    # Subset BAM to differential peaks only
    bedtools intersect -abam ${dedup_bam} -b ${differential_peaks_bed} > ${sample.id}_diffbind.bam
    samtools index ${sample.id}_diffbind.bam

    # Tumor-only Mutect2 with germline + PoN; --native-pair-hmm-threads 8 matches original 6_gatk-mutect2.sh
    gatk Mutect2 \\
        -R ${fasta} \\
        -I ${sample.id}_diffbind.bam \\
        --native-pair-hmm-threads 8 \\
        --germline-resource ${gnomad} \\
        --panel-of-normals ${pon} \\
        -O ${sample.id}.vcf

    # Atomize multi-allelic sites
    bcftools norm --atomize -o ${sample.id}_decomposed.vcf ${sample.id}.vcf

    # Drop alt-chromosome variants (chrUn_*, chr*_alt, etc.)
    awk '/^#/ || \$1 !~ /_/' ${sample.id}_decomposed.vcf > ${sample.id}_no-alt-chr.vcf

    # Filter C>T (matches original 6_gatk-mutect2.sh awk)
    awk 'BEGIN{OFS="\\t"} {
        if (\$0 ~ /^#/) { print; next }
        ref=\$4; alt=\$5
        if ((index(ref,"C") && index(alt,"T")) || (index(ref,"G") && index(alt,"A"))) print
    }' ${sample.id}_no-alt-chr.vcf > ${sample.id}_CT.vcf

    # Filter A>G
    awk 'BEGIN{OFS="\\t"} {
        if (\$0 ~ /^#/) { print; next }
        ref=\$4; alt=\$5
        if ((index(ref,"A") && index(alt,"G")) || (index(ref,"T") && index(alt,"C"))) print
    }' ${sample.id}_no-alt-chr.vcf > ${sample.id}_AG.vcf
    """
}
