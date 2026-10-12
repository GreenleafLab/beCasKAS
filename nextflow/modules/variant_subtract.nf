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

    # Keep the highest-FORMAT/AF record per CHROM:POS:REF:ALT. Atomizing a Mutect2 MNP
    # copies its whole FORMAT onto each base, so one replicate can carry the same SNV
    # twice (MNP piece + standalone call). Duplicates share a POS, so buffer per POS.
    dedup_max_af() {
        awk -F'\\t' '
            function flush(   i) { for (i = 1; i <= n; i++) print rec[keys[i]]; n = 0; delete rec; delete af }
            /^#/ { print; next }
            {
                loc = \$1 ":" \$2
                if (loc != cur) { flush(); cur = loc }
                k = \$4 ">" \$5
                nf = split(\$9, fmt, ":"); ai = 0
                for (i = 1; i <= nf; i++) if (fmt[i] == "AF") ai = i
                split(\$10, val, ":")
                a = (ai && val[ai] != "." && val[ai] != "") ? val[ai] + 0 : -1
                if (!(k in rec)) { keys[++n] = k; rec[k] = \$0; af[k] = a }
                else if (a > af[k]) { rec[k] = \$0; af[k] = a }
            }
            END { flush() }'
    }

    prep_and_merge() {
        local label="\$1"; shift
        local gz=() spec s v

        for spec in "\$@"; do
            s="\${spec%%:*}"
            v="\${spec#*:}"
            printf '%s\\n' "\$s" > "\$s.rename.txt"
            bcftools reheader -s "\$s.rename.txt" -o "\$s.renamed.vcf" "\$v"
            # Split MNPs and multi-allelic records per replicate, before merging: merging
            # first lets one replicate's MNP pieces pair with another's standalone calls
            # and leaves REF>ALT,* rows that downstream SNV filters drop. The split turns
            # each '*' (overlapping-allele placeholder) into its own record; drop those.
            bcftools norm --atomize -m -any "\$s.renamed.vcf" \\
                | bcftools view -e 'ALT="*"' \\
                | dedup_max_af | bgzip -c > "\$s.vcf.gz"
            rm -f "\$s.renamed.vcf"
            tabix -f -p vcf "\$s.vcf.gz"
            gz+=("\$s.vcf.gz")
        done

        if [[ "\${#gz[@]}" -gt 1 ]]; then
            # -m none: records join only on identical REF/ALT, so rows stay biallelic
            bcftools merge -m none -O z -o "\${label}_merged.vcf.gz" "\${gz[@]}"
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

    # already atomized and biallelic (prep_and_merge), so no norm step here
    awk '/^#/ || \$1 !~ /_/' \\
        \${prefix}_relcontrol_filt.vcf > \${prefix}_relcontrol_filt_no-alt-chr.vcf

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
