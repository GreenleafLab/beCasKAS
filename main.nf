#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

include { TRIM_ALIGN }                                                       from './nextflow/modules/trim_align'
include { PEAK_CALL; PEAK_CALL_BROAD; KAS_STATISTICS; NORMALIZE_BIGWIG; QC_METRICS } from './nextflow/modules/analyze'
include { DEDUP_CONSENSUS }                                                  from './nextflow/modules/picard_gencore_dedup'
include { SPLIT_BAM_PSEUDOREPS }                                             from './nextflow/modules/split_bam_pseudoreps'
include { BUILD_DIFFBIND_SAMPLESHEET; DIFFERENTIAL_PEAKS }                   from './nextflow/modules/differential_peaks'
include { CALL_VARIANTS }                                                    from './nextflow/modules/mutect2_call_variants'
include { VARIANT_SUBTRACT }                                                 from './nextflow/modules/variant_subtract'
include { PROTOSPACER }                                                      from './nextflow/modules/protospacer'
include { GUIDESCAN_SEED; GUIDESCAN_ENUMERATE; resolveGuidescanParams }      from './nextflow/modules/guidescan'

// Sharp non-differential arm (--sharp_nondiff) renamed processes.
include { GUIDESCAN_SEED   as GUIDESCAN_SEED_NONDIFF   }                       from './nextflow/modules/guidescan'
include { CALL_VARIANTS    as CALL_VARIANTS_NONDIFF    }                       from './nextflow/modules/mutect2_call_variants'
include { VARIANT_SUBTRACT as VARIANT_SUBTRACT_NONDIFF }                       from './nextflow/modules/variant_subtract'
include { PROTOSPACER      as PROTOSPACER_NONDIFF      }                       from './nextflow/modules/protospacer'

// Broad non-differential arm (--broad_nondiff) renamed processes.
include { BROAD_REGIONS; ARM_OVERLAP }                    from './nextflow/modules/broad_regions'
include { GUIDESCAN_SEED   as GUIDESCAN_SEED_BROAD   }    from './nextflow/modules/guidescan'
include { CALL_VARIANTS    as CALL_VARIANTS_BROAD    }    from './nextflow/modules/mutect2_call_variants'
include { VARIANT_SUBTRACT as VARIANT_SUBTRACT_BROAD }    from './nextflow/modules/variant_subtract'
include { PROTOSPACER      as PROTOSPACER_BROAD      }    from './nextflow/modules/protospacer'

// Stage One
workflow KAS_ANALYZER {

    take:
    ch_samples     // tuple(sample, fastq_1, fastq_2)
    refdir         // the --refdir directory as a file()
    fasta          // <genome>.fasta
    fasta_fai      // <genome>.fasta.fai
    assembly       // KAS-Analyzer assembly id (hg38, hg19)

    main:
    // Per-condition replicate counts
    ch_samples
        .filter { sample, r1, r2 -> sample.condition != 'background' }
        .map { sample, r1, r2 -> sample.condition }
        .reduce([:]) { acc, cond -> acc[cond] = (acc[cond] ?: 0) + 1; acc }
        .set { ch_replicate_counts }

    // Trim and align
    TRIM_ALIGN(ch_samples, refdir, assembly)
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
    PEAK_CALL(ch_bed_for_peaks, assembly)
    ch_broad_peaks = Channel.empty()
    if (params.broad_nondiff) {
        PEAK_CALL_BROAD(ch_bed_for_peaks, assembly)
        ch_broad_peaks = PEAK_CALL_BROAD.out.peaks
    }

    // Cohort-level RPM ratio computation
    ch_summary_noninput = TRIM_ALIGN.out.summary.filter { sample, summary -> sample.condition != 'background' }
    ch_summaries = ch_summary_noninput.collect { sample, summary -> summary }
    ch_summary_labels = ch_summary_noninput.collect { sample, summary -> sample.id }
    KAS_STATISTICS(ch_summaries, ch_summary_labels)

    // Normalization
    ch_bg.filter { sample, bg -> sample.condition != 'background' }
         .collect { sample, bg -> bg }
         .set { all_bg }
    NORMALIZE_BIGWIG(all_bg, KAS_STATISTICS.out.ratios, KAS_STATISTICS.out.labels, assembly)

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
        ch_diff_rows.toSortedList   { a, b -> a[0] <=> b[0] },
        ch_peaks_by_id.toSortedList { a, b -> a[0] <=> b[0] }
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
    broad_peaks  = ch_broad_peaks          // tuple(sample, epic2 bed); empty unless --broad_nondiff
    dedup        = ch_dedup                // tuple(sample, bam, bai)
}

// Stage Two
workflow DIFFERENTIAL {

    take:
    ch_case_samplesheets   // tuple(case_id, csv)
    ch_dedup               // tuple(sample, bam, bai)
    ch_broad_peaks         // tuple(sample, epic2 bed); empty unless --broad_nondiff
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

    ch_diff_bed = DIFFERENTIAL_PEAKS.out.bed.branch { case_id, bed ->
        with_peaks: bed.size() > 0
        no_peaks  : true
    }
    ch_diff_bed.no_peaks.subscribe { case_id, bed ->
        log.warn "DIFFERENTIAL: case '${case_id}' has no peaks passing the FDR/LFC " +
                 "filter -- skipping CALL_VARIANTS, VARIANT_SUBTRACT, both guidescan " +
                 "stages and PROTOSPACER"
    }

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
        .combine(ch_diff_bed.with_peaks, by: 0)
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

    // Background-subtracted per-case variant sets (see variant_subtract.nf).
    ch_subtract_jobs = CALL_VARIANTS.out.vcf_no_alt
        .map { case_id, sample, vcf ->
            tuple(case_id, sample.condition == 'control' ? 'control' : 'treatment', sample.id, vcf)
        }
        .groupTuple(by: 0)
        .map { case_id, groups, sample_ids, vcfs ->
            // Sort by sample id: groupTuple emits in arrival order and bcftools merge
            def ci = (0..<groups.size()).findAll { groups[it] == 'control'   }.sort { sample_ids[it] }
            def ti = (0..<groups.size()).findAll { groups[it] == 'treatment' }.sort { sample_ids[it] }
            tuple(case_id,
                  ci.collect { sample_ids[it] }, ci.collect { vcfs[it] },
                  ti.collect { sample_ids[it] }, ti.collect { vcfs[it] })
        }
        .branch { case_id, ctrl_ids, ctrl_vcfs, trt_ids, trt_vcfs ->
            complete: ctrl_ids && trt_ids
            partial : true
        }
    ch_subtract_jobs.partial.subscribe { case_id, ctrl_ids, ctrl_vcfs, trt_ids, trt_vcfs ->
        log.warn "VARIANT_SUBTRACT: case '${case_id}' has ${ctrl_ids.size()} control and " +
                 "${trt_ids.size()} treatment sample(s); subtraction needs at least one of " +
                 "each -- skipping"
    }
    VARIANT_SUBTRACT(ch_subtract_jobs.complete)

    // Treatment sample ids per case
    ch_tumors_by_case = ch_dedup
        .filter { sample, bam, bai -> sample.condition != 'control' }
        .map    { sample, bam, bai -> tuple(sample.condition, sample.id) }
        .groupTuple()
        .map    { cond, ids -> tuple(cond, ids.sort(false)) }

    // Two-stage gRNA binding-site prediction
    if (params.guidescan_grna || params.guidescan_params) {
        def guidescanOverrides = params.guidescan_params
            ? new groovy.json.JsonSlurper().parseText(file(params.guidescan_params, checkIfExists: true).text)
            : [:]

        ch_guidescan_jobs = ch_diff_bed.with_peaks
            .map    { case_id, bed -> tuple(case_id, bed, resolveGuidescanParams(guidescanOverrides, case_id)) }
            .filter { case_id, bed, gs -> gs != null }
        GUIDESCAN_SEED(ch_guidescan_jobs, fasta, fasta_fai)

        ch_enumerate_sites = Channel.empty()
        if (!params.skip_guidescan_enumerate) {
            ch_no_seed = GUIDESCAN_SEED.out.no_seed_bed.branch { case_id, bed, gs ->
                with_residue: bed.size() > 0
                empty       : true
            }
            ch_no_seed.empty.subscribe { case_id, bed, gs ->
                log.info "GUIDESCAN: case '${case_id}' has no unanchored peaks -- skipping enumerate"
            }
            GUIDESCAN_ENUMERATE(
                ch_no_seed.with_residue
                    .map { case_id, bed, gs -> tuple(case_id, bed, gs.findAll { k, v -> k != 'all_anchors' }) },
                fasta,
                fasta_fai
            )
            ch_enumerate_sites = GUIDESCAN_ENUMERATE.out.binding_sites
        }

        // Cohort-level, so it waits on both guidescan stages and VARIANT_SUBTRACT for every case
        if (!params.skip_protospacer) {
            ch_protospacer_in = GUIDESCAN_SEED.out.binding_sites
                .join(ch_diff_bed.with_peaks)
                .join(VARIANT_SUBTRACT.out.no_alt)
                .join(ch_tumors_by_case)
                .join(ch_guidescan_jobs.map { case_id, bed, gs -> tuple(case_id, gs.grna) })
                .toList()
                .map { rows ->
                    tuple(rows.collect { it[0] }, rows.collect { it[1] }, rows.collect { it[2] },
                          rows.collect { it[3] }, rows.collect { it[4] }, rows.collect { it[5] })
                }

            ch_enum_in = ch_enumerate_sites
                .toList()
                .map { rows -> tuple(rows.collect { it[0] }, rows.collect { it[1] }) }

            PROTOSPACER(
                ch_protospacer_in,
                ch_enum_in,
                fasta,
                fasta_fai,
                false
            )
        }

        // For highly transcribed loci (high endogenous ssDNA) that cannot clear DiffBind thresholds.
        if (params.sharp_nondiff) {
            ch_nondiff_jobs = DIFFERENTIAL_PEAKS.out.nondiff_peaks
                .filter { case_id, bed -> bed.size() > 0 }
                .map    { case_id, bed -> tuple(case_id, bed, resolveGuidescanParams(guidescanOverrides, case_id)) }
                .filter { case_id, bed, gs -> gs != null }
            GUIDESCAN_SEED_NONDIFF(ch_nondiff_jobs, fasta, fasta_fai)

            // Variant calling is restricted to peaks that got a scored site (seed)
            ch_nondiff_variant_jobs = ch_ctrl_by_case.mix(ch_trt_by_case)
                .combine(GUIDESCAN_SEED_NONDIFF.out.scored_peaks.filter { id, bed -> bed.size() > 0 }, by: 0)
                .map { case_id, sample, bam, bai, bed -> tuple(case_id, sample, bam, bai, bed) }
            CALL_VARIANTS_NONDIFF(
                ch_nondiff_variant_jobs,
                fasta, fasta_fai, fasta_dict, gnomad, gnomad_tbi, pon, pon_tbi
            )

            ch_nondiff_subtract = CALL_VARIANTS_NONDIFF.out.vcf_no_alt
                .map { case_id, sample, vcf ->
                    tuple(case_id, sample.condition == 'control' ? 'control' : 'treatment', sample.id, vcf)
                }
                .groupTuple(by: 0)
                .map { case_id, groups, sample_ids, vcfs ->
                    def ci = (0..<groups.size()).findAll { groups[it] == 'control'   }.sort { sample_ids[it] }
                    def ti = (0..<groups.size()).findAll { groups[it] == 'treatment' }.sort { sample_ids[it] }
                    tuple(case_id,
                          ci.collect { sample_ids[it] }, ci.collect { vcfs[it] },
                          ti.collect { sample_ids[it] }, ti.collect { vcfs[it] })
                }
                .filter { case_id, ctrl_ids, ctrl_vcfs, trt_ids, trt_vcfs -> ctrl_ids && trt_ids }
            VARIANT_SUBTRACT_NONDIFF(ch_nondiff_subtract)

            if (!params.skip_protospacer) {
                ch_nondiff_protospacer_in = GUIDESCAN_SEED_NONDIFF.out.binding_sites
                    .join(DIFFERENTIAL_PEAKS.out.nondiff_peaks)
                    .join(VARIANT_SUBTRACT_NONDIFF.out.no_alt)
                    .join(ch_tumors_by_case)
                    .join(ch_nondiff_jobs.map { case_id, bed, gs -> tuple(case_id, gs.grna) })
                    .toList()
                    .filter { rows -> rows.size() > 0 }
                    .map { rows ->
                        tuple(rows.collect { it[0] }, rows.collect { it[1] }, rows.collect { it[2] },
                              rows.collect { it[3] }, rows.collect { it[4] }, rows.collect { it[5] })
                    }
                // No enumerate in this arm, so the stage-2 channel is empty.
                PROTOSPACER_NONDIFF(
                    ch_nondiff_protospacer_in,
                    Channel.value(tuple([], [])),
                    fasta,
                    fasta_fai,
                    true
                )
            }
        }

        // Broad arm
        if (params.broad_nondiff) {
            ch_broad_ctrl_by_case = ch_broad_peaks
                .filter  { sample, peaks -> sample.condition == 'control' }
                .combine(ch_cases_only)
                .map     { sample, peaks, case_id -> tuple(case_id, peaks) }
            ch_broad_trt_by_case = ch_broad_peaks
                .filter  { sample, peaks -> sample.condition != 'control' }
                .map     { sample, peaks -> tuple(sample.condition, peaks) }
            // Sort by peakset name
            ch_broad_by_case = ch_broad_ctrl_by_case.mix(ch_broad_trt_by_case)
                .groupTuple()
                .map { case_id, ps -> tuple(case_id, ps.sort { it.name }) }

            ch_broad_min_overlap = ch_cases_only.map { case_id ->
                def o  = (diffbindParams[case_id] ?: [:])
                def mo = (o['dba.count']?.minOverlap != null) ? o['dba.count'].minOverlap
                                                              : params.min_overlap
                if (!(mo instanceof Integer) || mo < 1) {
                    error("case '${case_id}': dba.count.minOverlap=${mo} -- the broad arm needs " +
                          "an integer >= 1 (DiffBind's fractional minOverlap is not supported by the broad arm)")
                }
                tuple(case_id, mo as int)
            }

            BROAD_REGIONS(ch_broad_by_case.join(ch_broad_min_overlap))

            ch_broad_regions = BROAD_REGIONS.out.regions.branch { case_id, bed ->
                with_regions: bed.size() > 0
                empty       : true
            }
            ch_broad_regions.empty.subscribe { case_id, bed ->
                log.warn "BROAD_REGIONS: case '${case_id}' has no broad consensus regions -- skipping the broad arm"
            }
            ch_broad_jobs = ch_broad_regions.with_regions
                .map    { case_id, bed -> tuple(case_id, bed, resolveGuidescanParams(guidescanOverrides, case_id, true)) }
                .filter { case_id, bed, gs -> gs != null }
            GUIDESCAN_SEED_BROAD(ch_broad_jobs, fasta, fasta_fai)

            // Overlap with the other arms is reported here instead of subtracted from the broad regions.
            // join(remainder) emits the bare case_id, not (case_id, null), if the right side never emits.
            ch_sharp_scored = ch_cases_only
                .join(params.sharp_nondiff ? GUIDESCAN_SEED_NONDIFF.out.scored_peaks : Channel.empty(), remainder: true)
                .map { row -> row instanceof List ? tuple(row[0], row[1] ?: []) : tuple(row, []) }
            ARM_OVERLAP(
                DIFFERENTIAL_PEAKS.out.bed
                    .join(ch_sharp_scored)
                    .join(GUIDESCAN_SEED_BROAD.out.scored_peaks)
            )

            // variant calling is restricted to peaks with a scored site
            ch_broad_variant_jobs = ch_ctrl_by_case.mix(ch_trt_by_case)
                .combine(GUIDESCAN_SEED_BROAD.out.scored_peaks.filter { id, bed -> bed.size() > 0 }, by: 0)
                .map { case_id, sample, bam, bai, bed -> tuple(case_id, sample, bam, bai, bed) }
            CALL_VARIANTS_BROAD(
                ch_broad_variant_jobs,
                fasta, fasta_fai, fasta_dict, gnomad, gnomad_tbi, pon, pon_tbi
            )

            ch_broad_subtract = CALL_VARIANTS_BROAD.out.vcf_no_alt
                .map { case_id, sample, vcf ->
                    tuple(case_id, sample.condition == 'control' ? 'control' : 'treatment', sample.id, vcf)
                }
                .groupTuple(by: 0)
                .map { case_id, groups, sample_ids, vcfs ->
                    def ci = (0..<groups.size()).findAll { groups[it] == 'control'   }.sort { sample_ids[it] }
                    def ti = (0..<groups.size()).findAll { groups[it] == 'treatment' }.sort { sample_ids[it] }
                    tuple(case_id,
                          ci.collect { sample_ids[it] }, ci.collect { vcfs[it] },
                          ti.collect { sample_ids[it] }, ti.collect { vcfs[it] })
                }
                .filter { case_id, ctrl_ids, ctrl_vcfs, trt_ids, trt_vcfs -> ctrl_ids && trt_ids }
            VARIANT_SUBTRACT_BROAD(ch_broad_subtract)

            if (!params.skip_protospacer) {
                ch_broad_protospacer_in = GUIDESCAN_SEED_BROAD.out.binding_sites
                    .join(BROAD_REGIONS.out.regions)
                    .join(VARIANT_SUBTRACT_BROAD.out.no_alt)
                    .join(ch_tumors_by_case)
                    .join(ch_broad_jobs.map { case_id, bed, gs -> tuple(case_id, gs.grna) })
                    .toList()
                    .filter { rows -> rows.size() > 0 }
                    .map { rows ->
                        tuple(rows.collect { it[0] }, rows.collect { it[1] }, rows.collect { it[2] },
                              rows.collect { it[3] }, rows.collect { it[4] }, rows.collect { it[5] })
                    }
                // No enumerate in this arm, and the on-target site need not fall in a broad region.
                PROTOSPACER_BROAD(
                    ch_broad_protospacer_in,
                    Channel.value(tuple([], [])),
                    fasta,
                    fasta_fai,
                    true
                )
            }
        }
    }
}

// Workflow Entry
workflow {

    main:
    if (params.help) {
        log.info paramsHelp("nextflow run main.nf --samplesheet samples.csv --refdir /path/to/refs --outdir results",
                            parameters_schema: 'nextflow/configs/nextflow_schema.json')
        exit 0
    }

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
    def genome        = params.genome                         // reference file prefix, e.g. hg38, hg38_no_alt, hg19
    def assembly      = genome.split('_')[0]                  // assembly id for KAS-Analyzer and the GATK resources
    // KAS-Analyzer 2.0.0 ships blacklists/chrom sizes for these assemblies only
    def kas_assemblies = ['hg18', 'hg19', 'hg38', 'mm9', 'mm10', 'mm39', 'ce10', 'ce11', 'dm3', 'dm6', 'rn6', 'rn7', 'danRer10', 'danRer11']
    if (!(assembly in kas_assemblies)) {
        error("Invalid --genome '${genome}': '${assembly}' (the part before the first '_') must be a KAS-Analyzer assembly: ${kas_assemblies.join(', ')}.")
    }
    def fasta         = file("${refdir}/${genome}.fasta",                         checkIfExists: true)
    def fasta_fai     = file("${refdir}/${genome}.fasta.fai",                     checkIfExists: true)
    def fasta_dict    = file("${refdir}/${genome}.dict",                          checkIfExists: true)
    file("${refdir}/${genome}.1.bt2",                                             checkIfExists: true)
    def gnomad        = file("${refdir}/af-only-gnomad.${assembly}.vcf.gz",       checkIfExists: true)
    def gnomad_tbi    = file("${refdir}/af-only-gnomad.${assembly}.vcf.gz.tbi",   checkIfExists: true)
    def pon           = file("${refdir}/1000g_pon.${assembly}.vcf.gz",            checkIfExists: true)
    def pon_tbi       = file("${refdir}/1000g_pon.${assembly}.vcf.gz.tbi",        checkIfExists: true)

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
    KAS_ANALYZER(ch_samples, file(refdir), fasta, fasta_fai, assembly)

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
            KAS_ANALYZER.out.broad_peaks,
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
