#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

include { TRIM_ALIGN }                                                       from './nextflow/modules/trim_align'
include { PEAK_CALL; KAS_STATISTICS; NORMALIZE_BIGWIG; QC_METRICS }          from './nextflow/modules/analyze'
include { DEDUP_CONSENSUS }                                                  from './nextflow/modules/picard_gencore_dedup'
include { SPLIT_BAM_PSEUDOREPS }                                             from './nextflow/modules/split_bam_pseudoreps'
include { BUILD_DIFFBIND_SAMPLESHEET; DIFFERENTIAL_PEAKS }                   from './nextflow/modules/differential_peaks'
include { CALL_VARIANTS }                                                    from './nextflow/modules/mutect2_call_variants'

// Stage One
workflow KAS_ANALYZER {

    take:
    ch_samples     // tuple(sample, fastq_1, fastq_2)
    refdir         // the --refdir directory as a file()
    fasta          // hg38.fasta
    fasta_fai      // hg38.fasta.fai

    main:
    // Per-condition replicate counts
    ch_samples
        .filter { sample, r1, r2 -> sample.condition != 'background' }
        .map { sample, r1, r2 -> sample.condition }
        .reduce([:]) { acc, cond -> acc[cond] = (acc[cond] ?: 0) + 1; acc }
        .set { ch_replicate_counts }

    // Trim and align
    TRIM_ALIGN(ch_samples, refdir)
    ch_bam = TRIM_ALIGN.out.bam
    ch_bed = TRIM_ALIGN.out.bed
    ch_bg  = TRIM_ALIGN.out.bg

    // Peak calling
    ch_treatments = ch_bed
        .filter { sample, bed -> sample.condition != 'background' }
        .map    { sample, bed -> [sample.background_id, sample, bed] }
    ch_inputs = ch_bed
        .filter { sample, bed -> sample.condition == 'background' }
        .map    { sample, bed -> [sample.id, bed] }
    ch_bed_for_peaks = ch_treatments
        .combine(ch_inputs, by: 0)
        .map { ctrl_id, sample, sample_bed, ctrl_bed -> tuple(sample, sample_bed, ctrl_bed) }
    PEAK_CALL(ch_bed_for_peaks, 'hg38')

    // Cohort-level RPM ratio computation
    ch_summary_noninput = TRIM_ALIGN.out.summary.filter { sample, summary -> sample.condition != 'background' }
    ch_summaries = ch_summary_noninput.collect { sample, summary -> summary }
    ch_summary_labels = ch_summary_noninput.collect { sample, summary -> sample.id }
    KAS_STATISTICS(ch_summaries, ch_summary_labels)

    // Normalization
    ch_bg.filter { sample, bg -> sample.condition != 'background' }
         .collect { sample, bg -> bg }
         .set { all_bg }
    NORMALIZE_BIGWIG(all_bg, KAS_STATISTICS.out.ratios, KAS_STATISTICS.out.labels, 'hg38')

    // QC metrics - non-background samples only
    ch_bed_noninput = ch_bed.filter { sample, bed -> sample.condition != 'background' }
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

    // Picard/Gencode Deduplicate
    ch_bam_noninput = ch_bam.filter { sample, bam -> sample.condition != 'background' }
    DEDUP_CONSENSUS(ch_bam_noninput, fasta, fasta_fai)

    // Pseudo-replicate generation (split BAM into pseudo-reps if n=1)
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

    ch_branched.single
        .map { cond, samples, bams, bais -> tuple(samples[0], bams[0], bais[0]) }
        .set { ch_to_split }
    SPLIT_BAM_PSEUDOREPS(ch_to_split)

    ch_branched.multi
        .flatMap { cond, samples, bams, bais ->
            (0..<samples.size()).collect { i -> tuple(samples[i], bams[i], bais[i]) }
        }
        .set { ch_real_reps }

    // Setting up Differential Peaks samplesheet rows
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

    BUILD_DIFFBIND_SAMPLESHEET(
        ch_diff_rows.toList(),
        ch_peaks_by_id.toList()
    )

    // One (case_id, csv) per treatment case
    ch_case_samplesheets = BUILD_DIFFBIND_SAMPLESHEET.out.samplesheets
        .flatten()
        .map { ss ->
            def case_id = ss.name.replaceFirst(/^diffbind_samplesheet_/, '')
                                  .replaceFirst(/\.csv$/, '')
            tuple(case_id, ss)
        }

    emit:
    samplesheets = ch_case_samplesheets    // tuple(case_id, csv)
    dedup        = ch_dedup                // tuple(sample, bam, bai)
}

// Stage Two
workflow DIFFERENTIAL {

    take:
    ch_case_samplesheets   // tuple(case_id, csv)
    ch_dedup               // tuple(sample, bam, bai)
    fasta
    fasta_fai
    fasta_dict
    gnomad
    gnomad_tbi
    pon
    pon_tbi

    main:

    def diffbindParams = params.diffbind_params
        ? new groovy.json.JsonSlurper().parseText(file(params.diffbind_params, checkIfExists: true).text)
        : [:]
    ch_diff_jobs = ch_case_samplesheets.map { case_id, ss ->
        def o = (diffbindParams[case_id] ?: [:])
        def merged = [
            fdr     : (o.fdr      != null ? o.fdr      : params.fdr),
            lfc     : (o.lfc      != null ? o.lfc      : params.lfc),
            greylist: (o.greylist != null ? o.greylist : params.greylist)
        ]
        ['dba.count', 'dba.normalize', 'dba.analyze'].each { k ->
            if (o.containsKey(k)) merged[k] = o[k]
        }
        tuple(case_id, ss, groovy.json.JsonOutput.toJson(merged))
    }
    DIFFERENTIAL_PEAKS(ch_diff_jobs)

    // Per-case variant calling:
    ch_cases_only = ch_case_samplesheets.map { case_id, ss -> case_id }
    ch_ctrl_by_case = ch_dedup
        .filter { sample, bam, bai -> sample.condition == 'control' }
        .combine(ch_cases_only)
        .map { sample, bam, bai, case_id -> tuple(case_id, sample, bam, bai) }
    ch_trt_by_case = ch_dedup
        .filter { sample, bam, bai -> sample.condition != 'control' }
        .map { sample, bam, bai -> tuple(sample.condition, sample, bam, bai) }
    ch_variant_jobs = ch_ctrl_by_case.mix(ch_trt_by_case)
        .combine(DIFFERENTIAL_PEAKS.out.bed, by: 0)
        .map { case_id, sample, bam, bai, bed -> tuple(case_id, sample, bam, bai, bed) }

    CALL_VARIANTS(
        ch_variant_jobs,
        fasta,
        fasta_fai,
        fasta_dict,
        gnomad,
        gnomad_tbi,
        pon,
        pon_tbi
    )
}

// Workflow Entry
workflow {

    main:
    // Validate params
    validateParameters(parameters_schema: 'nextflow/configs/nextflow_schema.json')
    log.info paramsSummaryLog(workflow, parameters_schema: 'nextflow/configs/nextflow_schema.json')
    def valid_stages = ['all', 'kas_analyzer', 'diffbind']
    if (!(params.stage in valid_stages)) {
        error("Invalid --stage '${params.stage}'. Choose one of: ${valid_stages.join(', ')}.")
    }

    file(params.outdir).mkdirs()
    file("${params.outdir}/pipeline.running").text = ''

    // References from --refdir
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
        .map { sample_id, fastq_1, fastq_2, condition, replicate, background_id ->
            def sample = [
                id:            sample_id,
                condition:     condition,
                replicate:     replicate,
                background_id: background_id ?: ''
            ]
            tuple(sample, file(fastq_1), file(fastq_2))
        }
        .set { ch_samples }

    // Run Stage One
    KAS_ANALYZER(ch_samples, file(refdir), fasta, fasta_fai)

    // Stage Two (--stage diffbind)
    if (params.stage in ['all', 'diffbind']) {
        ch_case_samplesheets = params.diffbind_samplesheet_dir
            ? Channel.fromPath("${params.diffbind_samplesheet_dir}/diffbind_samplesheet_*.csv", checkIfExists: true)
                     .map { ss ->
                         def case_id = ss.name.replaceFirst(/^diffbind_samplesheet_/, '')
                                               .replaceFirst(/\.csv$/, '')
                         tuple(case_id, ss)
                     }
            : KAS_ANALYZER.out.samplesheets
        DIFFERENTIAL(
            ch_case_samplesheets,
            KAS_ANALYZER.out.dedup,
            fasta,
            fasta_fai,
            fasta_dict,
            gnomad,
            gnomad_tbi,
            pon,
            pon_tbi
        )
    }
}

workflow.onComplete {
    file("${params.outdir}/pipeline.running").delete()
    def marker = workflow.success ? 'pipeline.success' : 'pipeline.failure'
    file("${params.outdir}/${marker}").text = ''
}
