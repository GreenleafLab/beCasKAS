#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
})

# ---------- Helper functions (testable) ----------

#' Build a DiffBind-format sample sheet from the pipeline samplesheet.
build_diffbind_samplesheet <- function(pipeline_ss) {
    keep <- pipeline_ss$condition %in% c("control", "treatment")
    ss <- pipeline_ss[keep, , drop = FALSE]
    data.frame(
        SampleID   = ss$sample_id,
        Condition  = ss$condition,
        Replicate  = ss$replicate,
        bamReads   = ss$bam_path,
        Peaks      = ss$peaks_path,
        PeakCaller = "macs",
        stringsAsFactors = FALSE
    )
}

#' Apply Fold > lfc_cutoff AND FDR < fdr_cutoff filter.
filter_significant_peaks <- function(report, fdr_cutoff, lfc_cutoff) {
    keep <- report$Fold > lfc_cutoff & report$FDR < fdr_cutoff
    report[keep, , drop = FALSE]
}

#' Emit a human-readable record of the regime used (real reps vs pseudo).
write_method_metadata <- function(path, replicate_counts, n_pseudo, seed, fdr, lfc,
                                  is_pseudo_regime = FALSE) {
    lines <- c("# differential_peaks regime metadata", "")
    for (cond in names(replicate_counts)) {
        n <- replicate_counts[[cond]]
        if (is_pseudo_regime) {
            lines <- c(lines, sprintf(
                "%s: 1 sample -> %d pseudo-replicates (auto-split, seed=%d)",
                cond, n_pseudo, seed))
        } else {
            lines <- c(lines, sprintf(
                "%s: %d samples -> real biological replicates", cond, n))
        }
    }
    lines <- c(lines, "",
               sprintf("DiffBind ran with: minMembers=2, method=DESeq2, bGreylist=FALSE%s",
                       if (is_pseudo_regime) ", normalize=LIB (pseudo-rep fallback)" else ""),
               sprintf("Filter: Fold > %g AND FDR < %g", lfc, fdr),
               sprintf("Random seed: %d", seed))
    writeLines(lines, path)
}

resolve_diffbind_consts <- function(arglist) {
    lapply(arglist, function(v) {
        if (is.character(v) && length(v) == 1L && grepl("^DBA_", v)) {
            getExportedValue("DiffBind", v)
        } else {
            v
        }
    })
}

read_diffbind_params <- function(path) {
    if (is.null(path) || !nzchar(path)) return(list())
    jsonlite::fromJSON(path, simplifyVector = FALSE)
}

if (sys.nframe() == 0L && !interactive()) {

    option_list <- list(
        make_option("--samplesheet",   type = "character"),
        make_option("--output-dir",    type = "character", default = "."),
        make_option("--fdr",           type = "double",    default = 0.05),
        make_option("--lfc",           type = "double",    default = 1),
        make_option("--min-overlap",   type = "integer",   default = 2),
        make_option("--greylist",      type = "logical",   default = FALSE,
                    action = "store_true"),
        make_option("--seed",          type = "integer",   default = 42),
        make_option("--threads",       type = "integer",   default = 4),
        make_option("--n-pseudo-replicates", type = "integer", default = 3),
        make_option("--diffbind-params", type = "character", default = NULL)
    )
    args <- parse_args(OptionParser(option_list = option_list))
    set.seed(args$seed)

    suppressPackageStartupMessages({
        library(DiffBind)
        library(ggplot2)
    })

    # Read pipeline samplesheet
    pipeline_ss <- read.csv(args$samplesheet, stringsAsFactors = FALSE)
    dbp      <- read_diffbind_params(args$`diffbind-params`)
    fdr_cut  <- if (!is.null(dbp$fdr))      dbp$fdr      else args$fdr
    lfc_cut  <- if (!is.null(dbp$lfc))      dbp$lfc      else args$lfc
    greylist <- if (!is.null(dbp$greylist)) dbp$greylist else args$greylist

    # Compute per-condition rep counts for metadata
    rep_counts <- list(
        control   = sum(pipeline_ss$condition == "control"),
        treatment = sum(pipeline_ss$condition == "treatment")
    )

    # Build DiffBind sample sheet
    dba_ss <- build_diffbind_samplesheet(pipeline_ss)

    is_pseudo_regime <- any(grepl("_pseudo[0-9]+$", dba_ss$SampleID))

    # Run DiffBind
    samples <- dba(sampleSheet = dba_ss)
    count_args <- modifyList(
        list(minOverlap = args$`min-overlap`, bParallel = TRUE),
        resolve_diffbind_consts(dbp$`dba.count`))
    samples <- do.call(dba.count, c(list(samples), count_args))

    norm_args <- resolve_diffbind_consts(dbp$`dba.normalize`)
    if (is_pseudo_regime && is.null(norm_args$normalize)) {
        message("Pseudo-replicate regime detected; using library-size normalization (DBA_NORM_LIB).")
        norm_args$normalize <- DBA_NORM_LIB
    }
    if (length(norm_args)) {
        samples <- do.call(dba.normalize, c(list(samples), norm_args))
    }

    samples <- dba.contrast(samples, categories = DBA_CONDITION, minMembers = 2)
    analyze_args <- modifyList(
        list(bRetrieveAnalysis = FALSE, bParallel = TRUE, bGreylist = greylist),
        resolve_diffbind_consts(dbp$`dba.analyze`))
    samples <- do.call(dba.analyze, c(list(samples), analyze_args))
    dba.show(samples, bContrasts = TRUE)

    # Pull report
    report <- as.data.frame(dba.report(samples, contrast = 1, th = 1))

    # Color column for volcano
    report$color <- "NS"
    report$color[report$FDR < fdr_cut & report$Fold >  1] <- "Up"
    report$color[report$FDR < fdr_cut & report$Fold < 1]  <- "Down"

    # Volcano plot
    v <- ggplot(report, aes(x = Fold, y = -log10(FDR), color = color)) +
         geom_point(alpha = 0.2) +
         scale_color_manual(values = c("Up" = "#FCAE1E", "Down" = "grey70", "NS" = "grey70")) +
         theme_minimal() + xlim(-8, 8) + ylim(0, 12.5) +
         labs(x = "log2 Fold Change", y = "-log10 FDR") +
         theme(legend.position = "bottom")
    ggsave(file.path(args$`output-dir`, "volcano.pdf"), v, width = 7.5, height = 10)

    # Filtered differential_peaks.bed
    sig <- filter_significant_peaks(report,
                                    fdr_cutoff = fdr_cut,
                                    lfc_cutoff = lfc_cut)
    bed <- sig[, c("seqnames", "start", "end")]
    write.table(bed,
                file = file.path(args$`output-dir`, "differential_peaks.bed"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    # Full results TSV
    write.table(report,
                file = file.path(args$`output-dir`, "differential_peaks_full.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # Save DBA object for reproducibility
    saveRDS(samples, file = file.path(args$`output-dir`, "dba_object.rds"))

    # Regime metadata
    write_method_metadata(
        path             = file.path(args$`output-dir`, "differential_peaks_method.txt"),
        replicate_counts = rep_counts,
        n_pseudo         = args$`n-pseudo-replicates`,
        seed             = args$seed,
        fdr              = fdr_cut,
        lfc              = lfc_cut,
        is_pseudo_regime = is_pseudo_regime
    )

    cat("differential_peaks: wrote", nrow(bed), "significant peaks (FDR<", fdr_cut,
        ", Fold>", lfc_cut, ")\n")
}
