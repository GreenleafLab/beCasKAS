#!/usr/bin/env nextflow
// main.nf — beCasKAS pipeline entry point

nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-validation'

include { TRIM_ALIGN }                                                       from './nextflow/modules/trim_align'
include { PEAK_CALL; KAS_STATISTICS; NORMALIZE_BIGWIG; QC_METRICS }          from './nextflow/modules/analyze'
include { DEDUP_CONSENSUS }                                                  from './nextflow/modules/picard_gencore_dedup'
include { SPLIT_BAM_PSEUDOREPS }                                             from './nextflow/modules/split_bam_pseudoreps'
include { DIFFERENTIAL_PEAKS }                                               from './nextflow/modules/differential_peaks'
include { CALL_VARIANTS }                                                    from './nextflow/modules/mutect2_call_variants'

// Tool-version provenance: captures in-container tool/package versions that Nextflow's
// built-in tracing (.nextflow.log, execution_trace.txt, execution_report.html) cannot see.
// Output published to ${params.outdir}/pipeline_info/tool_versions.tsv.
process TOOL_VERSIONS {

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    path 'tool_versions.tsv'

    script:
    """
    set -euo pipefail
    {
        echo -e "tool\\tversion"
        echo -e "trim_galore\\t\$(trim_galore --version | grep -i 'version' | awk '{print \$NF}')"
        echo -e "bowtie2\\t\$(bowtie2 --version | head -1 | awk '{print \$NF}')"
        echo -e "samtools\\t\$(samtools --version | head -1 | awk '{print \$NF}')"
        echo -e "macs2\\t\$(macs2 --version 2>&1 | awk '{print \$NF}')"
        echo -e "picard\\t\$(picard MarkDuplicates --version 2>&1 | head -1)"
        echo -e "gencore\\t\$(gencore --version 2>&1 | head -1 | awk '{print \$NF}')"
        echo -e "gatk\\t\$(gatk --version | head -1 | awk '{print \$NF}')"
        echo -e "bcftools\\t\$(bcftools --version | head -1 | awk '{print \$NF}')"
        echo -e "bedtools\\t\$(bedtools --version | awk '{print \$NF}')"
        echo -e "deeptools\\t\$(multiBigwigSummary --version 2>&1 | awk '{print \$NF}')"
        echo -e "R\\t\$(R --version | head -1 | awk '{print \$3}')"
        echo -e "DiffBind\\t\$(Rscript -e 'cat(as.character(packageVersion(\\"DiffBind\\")))')"
        echo -e "DESeq2\\t\$(Rscript -e 'cat(as.character(packageVersion(\\"DESeq2\\")))')"
        echo -e "KAS-Analyzer\\t1.0"
    } > tool_versions.tsv
    """
}

process BUILD_DIFFBIND_SAMPLESHEET {

    input:
    val rows         // list of [sample_id, condition, replicate, bam_path]
    val peaks_map    // list of [sample_id, peaks_path] (real samples only; pseudo-reps reuse parent peaks)

    output:
    path 'diffbind_samplesheet.csv'

    script:
    """
    python3 - <<'PYEOF' > diffbind_samplesheet.csv
import csv, sys

rows = ${groovy.json.JsonOutput.toJson(rows)}
peaks_map = dict([(p[0], p[1]) for p in ${groovy.json.JsonOutput.toJson(peaks_map)}])

writer = csv.writer(sys.stdout)
writer.writerow(['sample_id', 'condition', 'replicate', 'bam_path', 'peaks_path'])
for r in rows:
    sid, cond, rep, bam = r
    # Pseudo-rep IDs look like "X_pseudo1"; their peaks are the parent's
    parent = sid.split('_pseudo')[0]
    peaks = peaks_map.get(sid) or peaks_map.get(parent) or ''
    writer.writerow([sid, cond, rep, bam, peaks])
PYEOF
    """
}

workflow {

    // 1. Validate params
    validateParameters(parameters_schema: 'nextflow/configs/nextflow_schema.json')
    log.info paramsSummaryLog(workflow, parameters_schema: 'nextflow/configs/nextflow_schema.json')
    TOOL_VERSIONS()

    // 2. Resolve references from --reference_dir (user-staged directory)
    if (!params.reference_dir) {
        error("Missing required --reference_dir. See README for required layout.")
    }
    def ref_dir       = params.reference_dir
    def fasta         = file("${ref_dir}/hg38.fasta",                       checkIfExists: true)
    def fasta_fai     = file("${ref_dir}/hg38.fasta.fai",                   checkIfExists: true)
    def fasta_dict    = file("${ref_dir}/hg38.dict",                        checkIfExists: true)
    def bowtie2_index = file("${ref_dir}/hg38_Bowtie2Index",                checkIfExists: true)
    def gnomad        = file("${ref_dir}/af-only-gnomad.hg38.vcf.gz",       checkIfExists: true)
    def gnomad_tbi    = file("${ref_dir}/af-only-gnomad.hg38.vcf.gz.tbi",   checkIfExists: true)
    def pon           = file("${ref_dir}/1000g_pon.hg38.vcf.gz",            checkIfExists: true)
    def pon_tbi       = file("${ref_dir}/1000g_pon.hg38.vcf.gz.tbi",        checkIfExists: true)

    // 3. Parse samplesheet → channels
    Channel
        .fromSamplesheet('samplesheet', skip_duplicate_check: false)
        .map { row ->
            def sample = [
                id:         row.sample_id,
                condition:  row.condition,
                replicate:  row.replicate,
                control_id: row.control_id ?: ''
            ]
            tuple(sample, file(row.fastq_1), file(row.fastq_2))
        }
        .set { ch_samples }

    // Per-condition replicate counts (drives pseudo-rep dispatch later)
    ch_samples
        .filter { sample, r1, r2 -> sample.condition != 'input' }
        .map { sample, r1, r2 -> sample.condition }
        .reduce([:]) { acc, cond -> acc[cond] = (acc[cond] ?: 0) + 1; acc }
        .set { ch_replicate_counts }

    // 4. trim_align — every sample
    TRIM_ALIGN(ch_samples, file(ref_dir))

    // 5. Build channels for downstream modules

    // Per-sample BAM/BED/BG channels (joined by sample)
    ch_bam = TRIM_ALIGN.out.bam     // tuple(sample, sorted.bam)
    ch_bed = TRIM_ALIGN.out.bed     // tuple(sample, .bed)
    ch_bg  = TRIM_ALIGN.out.bg      // tuple(sample, .bg)

    // 6. Peak calling — non-input samples paired with their control_id BED
    ch_bed_by_id = ch_bed.map { sample, bed -> [sample.id, sample, bed] }
    ch_bed_for_peaks = ch_bed_by_id
        .filter { id, sample, bed -> sample.condition != 'input' }
        .map { id, sample, bed -> [sample.control_id, sample, bed] }
        .combine(ch_bed_by_id.map { id, sample, bed -> [id, bed] }, by: 0)
        .map { ctrl_id, sample, sample_bed, ctrl_bed -> tuple(sample, sample_bed, ctrl_bed) }

    PEAK_CALL(ch_bed_for_peaks, 'hg38')

    // 7a. Cohort-level RPM ratio computation (verbatim awk from 3_analyze.sh:74-96)
    // Collect per-sample Summary/ dirs from TRIM_ALIGN, run KAS-Analyzer statistics, awk → ratios.txt
    ch_summaries = TRIM_ALIGN.out.summary.collect { sample, summary -> summary }
    ch_summary_labels = TRIM_ALIGN.out.summary.collect { sample, summary -> sample.id }

    KAS_STATISTICS(ch_summaries, ch_summary_labels)

    // 7b. Normalization
    ch_bg.collect { sample, bg -> bg }.set { all_bg }
    NORMALIZE_BIGWIG(all_bg, KAS_STATISTICS.out.ratios, KAS_STATISTICS.out.labels, 'hg38')

    // 8. QC metrics — non-input only (matches original)
    ch_bed_noninput = ch_bed.filter { sample, bed -> sample.condition != 'input' }
    ch_peaks_for_qc = PEAK_CALL.out.peaks_abridged   // headerless BED for FRiP (bedtools-based)
    ch_bigwigs_for_qc = NORMALIZE_BIGWIG.out.bigwigs.flatten().filter {
        // exclude input samples' bigwigs from correlation
        true   // TODO: filter by sample id once channels are joined
    }

    QC_METRICS(
        ch_bed_noninput.collect { sample, bed -> bed },
        ch_peaks_for_qc.collect { sample, peaks -> peaks },
        ch_bigwigs_for_qc.collect(),
        ch_bed_noninput.collect { sample, bed -> sample.id }
    )

    // 9. Dedup — non-input only
    ch_bam_noninput = ch_bam.filter { sample, bam -> sample.condition != 'input' }
    DEDUP_CONSENSUS(ch_bam_noninput, fasta, fasta_fai)

    // 10. Pseudo-rep dispatch: per condition, if n=1, split that sample's BAM into 3 pseudo-reps
    ch_dedup = DEDUP_CONSENSUS.out.bam.join(DEDUP_CONSENSUS.out.bai)
        // tuple(sample, dedup_bam, dedup_bai)

    // Build a map condition → list of (sample, bam, bai) entries
    ch_dedup_by_cond = ch_dedup.map { sample, bam, bai -> [sample.condition, sample, bam, bai] }

    // Group by condition; if size==1 → split; else pass through as-is
    ch_dedup_by_cond
        .groupTuple()    // [cond, [sample...], [bam...], [bai...]]
        .branch {
            single:    it[1].size() == 1
            multi:     it[1].size() >= 2
            none:      true
        }
        .set { ch_branched }

    // Single-rep conditions → SPLIT_BAM_PSEUDOREPS
    ch_branched.single
        .map { cond, metas, bams, bais ->
            tuple(metas[0], bams[0], bais[0])
        }
        .set { ch_to_split }

    SPLIT_BAM_PSEUDOREPS(ch_to_split)

    // Multi-rep conditions → flatten back to per-sample tuples
    ch_branched.multi
        .flatMap { cond, metas, bams, bais ->
            (0..<metas.size()).collect { i -> tuple(metas[i], bams[i], bais[i]) }
        }
        .set { ch_real_reps }

    // Build the unified samplesheet for DIFFERENTIAL_PEAKS:
    //   each row has sample_id, condition, replicate, bam_path, peaks_path
    //   pseudo-reps get synthesized rows (sample_id = original + "_pseudoN")
    ch_real_reps_for_diff = ch_real_reps.map { sample, bam, bai ->
        [sample.id, sample.condition, sample.replicate, bam.toString()]
    }

    SPLIT_BAM_PSEUDOREPS.out.pseudoreps
        .flatMap { sample, bams, bais ->
            (0..<bams.size()).collect { i ->
                [sample.id + "_pseudo${i+1}", sample.condition, i+1, bams[i].toString()]
            }
        }
        .set { ch_pseudo_for_diff }

    // Combine real + pseudo, then join with peaks
    ch_diff_rows = ch_real_reps_for_diff.mix(ch_pseudo_for_diff)

    PEAK_CALL.out.peaks
        .map { sample, peaks -> [sample.id, peaks.toString()] }
        .set { ch_peaks_by_id }

    // Build samplesheet TSV (one row per BAM, joined to peaks)
    BUILD_DIFFBIND_SAMPLESHEET(
        ch_diff_rows.collect(),
        ch_peaks_by_id.collect()
    )
        .set { ch_diff_samplesheet }

    DIFFERENTIAL_PEAKS(ch_diff_samplesheet)

    // 11. call_variants — per non-input sample, intersected with cohort differential_peaks.bed
    ch_dedup_for_variants = DEDUP_CONSENSUS.out.bam.join(DEDUP_CONSENSUS.out.bai)

    CALL_VARIANTS(
        ch_dedup_for_variants,
        DIFFERENTIAL_PEAKS.out.bed,
        fasta,
        fasta_fai,
        fasta_dict,
        gnomad,
        gnomad_tbi,
        pon,
        pon_tbi
    )
}
