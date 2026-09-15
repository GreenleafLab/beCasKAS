// Background-subtracted per-case variant set, from CALL_VARIANTS' per-sample VCFs.

process VARIANT_SUBTRACT {

    tag "${case_id}"

    input:
    tuple val(case_id), val(control_ids), path(control_vcfs), val(treatment_ids), path(treatment_vcfs)

    output:
    tuple val(case_id), path("${case_id}_vs_control_relcontrol_filt_no-alt-chr.vcf"), emit: no_alt
    tuple val(case_id), path("${case_id}_vs_control_relcontrol_filt_CT.vcf"),         emit: vcf_ct
    tuple val(case_id), path("${case_id}_vs_control_relcontrol_filt_AG.vcf"),         emit: vcf_ag

    script:
    def control_spec   = [control_ids,   control_vcfs  ].transpose().collect { id, v -> "${id}:${v}" }.join(' ')
    def treatment_spec = [treatment_ids, treatment_vcfs].transpose().collect { id, v -> "${id}:${v}" }.join(' ')
    """
    set -euo pipefail

    prep_and_merge() {
        local label="\$1"; shift
        local gz=() spec s v

        for spec in "\$@"; do
            s="\${spec%%:*}"
            v="\${spec#*:}"
            printf '%s\\n' "\$s" > "\$s.rename.txt"
            bcftools reheader -s "\$s.rename.txt" -o "\$s.renamed.vcf" "\$v"
            bgzip -f "\$s.renamed.vcf"
            mv "\$s.renamed.vcf.gz" "\$s.vcf.gz"
            tabix -f -p vcf "\$s.vcf.gz"
            gz+=("\$s.vcf.gz")
        done

        if [[ "\${#gz[@]}" -gt 1 ]]; then
            bcftools merge -O z -o "\${label}_merged.vcf.gz" "\${gz[@]}"
        else
            cp "\${gz[0]}" "\${label}_merged.vcf.gz"
        fi
        tabix -f -p vcf "\${label}_merged.vcf.gz"
    }

    prep_and_merge control ${control_spec}
    prep_and_merge treatment ${treatment_spec}

    prefix=${case_id}_vs_control

    bcftools isec -C -w1 treatment_merged.vcf.gz control_merged.vcf.gz \\
        -O v -o \${prefix}_relcontrol.vcf

    bcftools view -i 'FORMAT/AD[*:1] >= 2 && FORMAT/DP[*] >= 5' \\
        \${prefix}_relcontrol.vcf -o \${prefix}_relcontrol_filt.vcf

    bcftools norm --atomize \\
        -o \${prefix}_relcontrol_filt_decomposed.vcf \${prefix}_relcontrol_filt.vcf

    awk '/^#/ || \$1 !~ /_/' \\
        \${prefix}_relcontrol_filt_decomposed.vcf > \${prefix}_relcontrol_filt_no-alt-chr.vcf

    awk '{
        if (\$0 ~ /^#/) { print; next }
        ref = \$4; alt = \$5
        if ((index(ref, "C") && index(alt, "T")) || (index(ref, "G") && index(alt, "A"))) print
    }' \${prefix}_relcontrol_filt_no-alt-chr.vcf > \${prefix}_relcontrol_filt_CT.vcf

    awk '{
        if (\$0 ~ /^#/) { print; next }
        ref = \$4; alt = \$5
        if ((index(ref, "A") && index(alt, "G")) || (index(ref, "T") && index(alt, "C"))) print
    }' \${prefix}_relcontrol_filt_no-alt-chr.vcf > \${prefix}_relcontrol_filt_AG.vcf

    # `|| true` because grep -c exits 1 on zero matches, a legitimate result here.
    echo "=== ${case_id} vs control (merged replicates + subtract): \\
\$(zgrep -vc '^#' treatment_merged.vcf.gz || true) raw treatment, \\
\$(grep -vc '^#' \${prefix}_relcontrol.vcf || true) post-subtract, \\
\$(grep -vc '^#' \${prefix}_relcontrol_filt.vcf || true) post-filter, \\
\$(grep -vc '^#' \${prefix}_relcontrol_filt_CT.vcf || true) CT, \\
\$(grep -vc '^#' \${prefix}_relcontrol_filt_AG.vcf || true) AG ==="
    """
}
