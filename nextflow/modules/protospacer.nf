// Variant position within the gRNA protospacer, per differential peak. See the
// README section "Protospacer Position Analysis".
//
// COHORT-LEVEL, unlike every other process in this stage: same-locus groups share a
// y-axis and the funnel figure is one combined plot, so it needs every case at once
// rather than one task per case. Runs once per STRANDED mode.

process PROTOSPACER {

    tag "${stranded ? 'stranded' : 'aggregated'}"

    input:
    tuple val(case_ids), path(seed_tsvs, stageAs: 'seed_??.tsv'), path(diff_beds, stageAs: 'diff_??.bed'), path(sub_vcfs), val(tumor_lists), val(grnas)
    tuple val(enum_ids), path(enum_tsvs, stageAs: 'enum_??.tsv')
    val  stranded
    path fasta
    path fasta_fai

    output:
    path "out/**", emit: results

    script:
    // Largest same-spacer set shares a y-axis; cross-locus cases get their own.
    def by_grna = [case_ids, grnas].transpose().groupBy { id, g -> g }
    def shared  = (by_grna.sort { a, b -> b.value.size() <=> a.value.size() ?: a.key <=> b.key }
                          .values().find { true } ?: []).collect { id, g -> id }

    def ordered = case_ids.sort(false)
    def ramp    = ['#D97827', '#1565C0', '#BF2025', '#F4B97D', '#0B6623', '#6A4C93', '#FCAE1E', '#00838F']
    def palette = params.protospacer_palette
        ? new groovy.json.JsonSlurper().parseText(file(params.protospacer_palette, checkIfExists: true).text)
        : ordered.withIndex().collectEntries { id, i -> [id, ramp[i % ramp.size()]] }

    def n = ordered.size()
    def cfg = [
        groups          : ordered,
        tumor_samples   : [case_ids, tumor_lists].transpose().collectEntries { id, t -> [id, t] },
        seed_files      : ordered.collectEntries { [it, "${it}_binding_sites_seed.tsv"] },
        enumerate_files : ordered.collectEntries { [it, "${it}_binding_sites_enumerate.tsv"] },
        peakbed_files   : ordered.collectEntries { [it, "${it}_differential_peaks.bed"] },
        palette         : palette,
        shared_groups   : shared,
        panel_ncol      : n <= 3 ? n  : 2,
        panel_w         : n <= 3 ? 16 : 15,
        panel_h         : n <= 3 ? 5.2 : 9
    ]

    def seed_spec = [case_ids, seed_tsvs].transpose().collect { id, f -> "${id}:${f}" }.join(' ')
    def bed_spec  = [case_ids, diff_beds].transpose().collect { id, f -> "${id}:${f}" }.join(' ')
    def vcf_spec  = [case_ids, sub_vcfs ].transpose().collect { id, f -> "${id}:${f}" }.join(' ')
    def enum_spec = enum_ids ? [enum_ids, enum_tsvs].transpose().collect { id, f -> "${id}:${f}" }.join(' ') : ''
    """
    set -euo pipefail

    Rscript -e 'if (!requireNamespace("patchwork", quietly=TRUE)) stop(
      "PROTOSPACER needs r-patchwork", call.=FALSE)'

    cat > protospacer_config.json <<'JSON'
    ${groovy.json.JsonOutput.toJson(cfg)}
    JSON

    mkdir -p sites diff vcf out

    for spec in ${seed_spec}; do
        cp "\${spec#*:}" "sites/\${spec%%:*}_binding_sites_seed.tsv"
    done
    for spec in ${bed_spec}; do
        cp "\${spec#*:}" "diff/\${spec%%:*}_differential_peaks.bed"
    done

    # protospacer.R reads <vcf_dir>/<group>/<group>_vs_control_...vcf
    for spec in ${vcf_spec}; do
        c="\${spec%%:*}"
        mkdir -p "vcf/\$c"
        cp "\${spec#*:}" "vcf/\$c/\${c}_vs_control_relcontrol_filt_no-alt-chr.vcf"
    done
    for spec in ${enum_spec}; do
        cp "\${spec#*:}" "sites/\${spec%%:*}_binding_sites_enumerate.tsv"
    done

    ${projectDir}/nextflow/bin/protospacer.R \\
        --config protospacer_config.json \\
        --stranded ${stranded ? 'TRUE' : 'FALSE'} \\
        --seed-dir sites \\
        --vcf-dir vcf \\
        --diff-dir diff \\
        --out-dir out \\
        --genome ${fasta} \\
        --samtools \$(command -v samtools)

    # bgzip + tabix every VCF the stranded pass wrote
    for v in out/vcf/*.vcf; do
        [ -e "\$v" ] || break
        bgzip -c "\$v" > "\$v.gz"
        tabix -f -p vcf "\$v.gz"
    done
    """
}
