#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

include { TRIM_ALIGN }                                                       from './nextflow/modules/trim_align'
include { PEAK_CALL; KAS_STATISTICS; NORMALIZE_BIGWIG; QC_METRICS }          from './nextflow/modules/analyze'
include { DEDUP_CONSENSUS }                                                  from './nextflow/modules/picard_gencore_dedup'
include { SPLIT_BAM_PSEUDOREPS }                                             from './nextflow/modules/split_bam_pseudoreps'
include { DIFFERENTIAL_PEAKS }                                               from './nextflow/modules/differential_peaks'
include { CALL_VARIANTS }                                                    from './nextflow/modules/mutect2_call_variants'

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

    // Validate params
    validateParameters(parameters_schema: 'nextflow/configs/nextflow_schema.json')
    log.info paramsSummaryLog(workflow, parameters_schema: 'nextflow/configs/nextflow_schema.json')

    // Get references from --refdir
    if (!params.refdir) {
        error("Missing required --refdir. See README for required layout.")
    }
    def refdir        = params.refdir
    def fasta         = file("${refdir}/hg38.fasta",                       checkIfExists: true)
    def fasta_fai     = file("${refdir}/hg38.fasta.fai",                   checkIfExists: true)
    def fasta_dict    = file("${refdir}/hg38.dict",                        checkIfExists: true)
    file("${refdir}/hg38.1.bt2",                                           checkIfExists: true)
    def gnomad        = file("${refdir}/af-only-gnomad.hg38.vcf.gz",       checkIfExists: true)
    def gnomad_tbi    = file("${refdir}/af-only-gnomad.hg38.vcf.gz.tbi",   checkIfExists: true)
    def pon           = file("${refdir}/1000g_pon.hg38.vcf.gz",            checkIfExists: true)
    def pon_tbi       = file("${refdir}/1000g_pon.hg38.vcf.gz.tbi",        checkIfExists: true)

    // Parse samplesheet
    Channel
        .fromSamplesheet('samplesheet', parameters_schema: 'nextflow/configs/nextflow_schema.json', skip_duplicate_check: false)
        .map { sample_id, fastq_1, fastq_2, condition, replicate, control_id ->
            def sample = [
                id:         sample_id,
                condition:  condition,
                replicate:  replicate,
                control_id: control_id ?: ''
            ]
            tuple(sample, file(fastq_1), file(fastq_2))
        }
        .set { ch_samples }

    // Per-condition replicate counts
    ch_samples
        .filter { sample, r1, r2 -> sample.condition != 'input' }
        .map { sample, r1, r2 -> sample.condition }
        .reduce([:]) { acc, cond -> acc[cond] = (acc[cond] ?: 0) + 1; acc }
        .set { ch_replicate_counts }

    // Trim and align
    TRIM_ALIGN(ch_samples, file(refdir))
    ch_bam = TRIM_ALIGN.out.bam
    ch_bed = TRIM_ALIGN.out.bed
    ch_bg  = TRIM_ALIGN.out.bg

    // Peak calling
    ch_bed_by_id = ch_bed.map { sample, bed -> [sample.id, sample, bed] }
    ch_bed_for_peaks = ch_bed_by_id
        .filter { id, sample, bed -> sample.condition != 'input' }
        .map { id, sample, bed -> [sample.control_id, sample, bed] }
        .combine(ch_bed_by_id.map { id, sample, bed -> [id, bed] }, by: 0)
        .map { ctrl_id, sample, sample_bed, ctrl_bed -> tuple(sample, sample_bed, ctrl_bed) }
    PEAK_CALL(ch_bed_for_peaks, 'hg38')

    // Cohort-level RPM ratio computation
    ch_summaries = TRIM_ALIGN.out.summary.collect { sample, summary -> summary }
    ch_summary_labels = TRIM_ALIGN.out.summary.collect { sample, summary -> sample.id }
    KAS_STATISTICS(ch_summaries, ch_summary_labels)

    // Normalization
    ch_bg.collect { sample, bg -> bg }.set { all_bg }
    NORMALIZE_BIGWIG(all_bg, KAS_STATISTICS.out.ratios, KAS_STATISTICS.out.labels, 'hg38')

    // QC metrics -- non-input only (matches original)
    ch_bed_noninput = ch_bed.filter { sample, bed -> sample.condition != 'input' }
    ch_peaks_for_qc = PEAK_CALL.out.peaks_abridged
    ch_bigwigs_for_qc = NORMALIZE_BIGWIG.out.bigwigs.flatten().filter {
        true   // TODO: filter by sample id once channels are joined
    }

    QC_METRICS(
        ch_bed_noninput.collect { sample, bed -> bed },
        ch_peaks_for_qc.collect { sample, peaks -> peaks },
        ch_bigwigs_for_qc.collect(),
        ch_bed_noninput.collect { sample, bed -> sample.id }
    )

    // Deduplicate
    ch_bam_noninput = ch_bam.filter { sample, bam -> sample.condition != 'input' }
    DEDUP_CONSENSUS(ch_bam_noninput, fasta, fasta_fai)

    // Pseudo-replicate generation (split BAM into 3 pseudo-reps if n=1)
    ch_dedup = DEDUP_CONSENSUS.out.bam.join(DEDUP_CONSENSUS.out.bai)
    ch_dedup_by_cond = ch_dedup.map { sample, bam, bai -> [sample.condition, sample, bam, bai] }
    ch_dedup_by_cond
        .groupTuple()
        .branch {
            single:    it[1].size() == 1
            multi:     it[1].size() >= 2
            none:      true
        }
        .set { ch_branched }

    // Single-rep conditions
    ch_branched.single
        .map { cond, samples, bams, bais ->
            tuple(samples[0], bams[0], bais[0])
        }
        .set { ch_to_split }
    SPLIT_BAM_PSEUDOREPS(ch_to_split)

    // Multi-rep conditions
    ch_branched.multi
        .flatMap { cond, samples, bams, bais ->
            (0..<samples.size()).collect { i -> tuple(samples[i], bams[i], bais[i]) }
        }
        .set { ch_real_reps }

    // Build the unified samplesheet for DIFFERENTIAL_PEAKS:
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

    // Call Variants
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
