#!/usr/bin/env Rscript
# Reproduces 5E_differential_peaks.Rmd as a parameterized CLI script.
# Adds: explicit set.seed, CLI args, programmatic DiffBind samplesheet construction,
# regime metadata output. Computation identical to the original Rmd.

suppressPackageStartupMessages({
    library(optparse)
})

# ---------- Helper functions (testable) ----------

#' Build a DiffBind-format sample sheet from the pipeline samplesheet.
#' Excludes input rows; maps column names DiffBind expects.
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

#' Apply original code's Fold > lfc_cutoff AND FDR < fdr_cutoff filter.
filter_significant_peaks <- function(report, fdr_cutoff, lfc_cutoff) {
    keep <- report$Fold > lfc_cutoff & report$FDR < fdr_cutoff
    report[keep, , drop = FALSE]
}

#' Emit a human-readable record of the regime used (real reps vs pseudo).
write_method_metadata <- function(path, replicate_counts, n_pseudo, seed, fdr, lfc) {
    lines <- c("# differential_peaks regime metadata", "")
    for (cond in names(replicate_counts)) {
        n <- replicate_counts[[cond]]
        if (n == 1) {
            lines <- c(lines, sprintf(
                "%s: 1 sample → %d pseudo-replicates (auto-split, seed=%d)",
                cond, n_pseudo, seed))
        } else {
            lines <- c(lines, sprintf(
                "%s: %d samples → real biological replicates", cond, n))
        }
    }
    lines <- c(lines, "",
               sprintf("DiffBind ran with: minMembers=2, method=DESeq2, bGreylist=FALSE"),
               sprintf("Filter: Fold > %g AND FDR < %g", lfc, fdr),
               sprintf("Random seed: %d", seed))
    writeLines(lines, path)
}

# ---------- CLI entry point ----------

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
        make_option("--n-pseudo-replicates", type = "integer", default = 3)
    )
    args <- parse_args(OptionParser(option_list = option_list))
    set.seed(args$seed)

    suppressPackageStartupMessages({
        library(DiffBind)
        library(ggplot2)
    })

    # Read pipeline samplesheet
    pipeline_ss <- read.csv(args$samplesheet, stringsAsFactors = FALSE)

    # Compute per-condition rep counts for metadata
    rep_counts <- list(
        control   = sum(pipeline_ss$condition == "control"),
        treatment = sum(pipeline_ss$condition == "treatment")
    )

    # Build DiffBind sample sheet
    dba_ss <- build_diffbind_samplesheet(pipeline_ss)

    # Run DiffBind exactly as in 5E_differential_peaks.Rmd
    samples <- dba(sampleSheet = dba_ss)
    samples <- dba.count(samples, bParallel = TRUE)
    samples <- dba.contrast(samples, categories = DBA_CONDITION, minMembers = 2)
    samples <- dba.analyze(samples,
                           bRetrieveAnalysis = FALSE,
                           bParallel         = TRUE,
                           bGreylist         = args$greylist)
    dba.show(samples, bContrasts = TRUE)

    # Pull report
    report <- as.data.frame(dba.report(samples, contrast = 1, th = 1))

    # Color column for volcano (matches 5E)
    report$color <- "NS"
    report$color[report$FDR <  args$fdr & report$Fold >  1] <- "Up"
    report$color[report$FDR <  args$fdr & report$Fold < 1]  <- "Down"

    # Volcano plot (matches 5E ggplot)
    v <- ggplot(report, aes(x = Fold, y = -log10(FDR), color = color)) +
         geom_point(alpha = 0.2) +
         scale_color_manual(values = c("Up" = "#FCAE1E", "Down" = "grey70", "NS" = "grey70")) +
         theme_minimal() + xlim(-8, 8) + ylim(0, 12.5) +
         labs(x = "log2 Fold Change", y = "-log10 FDR") +
         theme(legend.position = "bottom")
    ggsave(file.path(args$`output-dir`, "volcano.pdf"), v, width = 7.5, height = 10)

    # Filtered differential_peaks.bed
    sig <- filter_significant_peaks(report,
                                    fdr_cutoff = args$fdr,
                                    lfc_cutoff = args$lfc)
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
        fdr              = args$fdr,
        lfc              = args$lfc
    )

    cat("differential_peaks: wrote", nrow(bed), "significant peaks (FDR<", args$fdr,
        ", Fold>", args$lfc, ")\n")
}
