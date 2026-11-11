#!/usr/bin/env Rscript
# Bins background-subtracted SNVs by strand-correct protospacer position (1 = PAM-distal)
# and writes the per-condition tables, figures and protospacer-edit VCFs.

suppressMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# Helvetica = Arial metrics in the PDFs (Arial isn't in the container); PNGs fall
# back to DejaVu Sans, so use the PDFs for figures.
PLOT_FONT <- "Helvetica"
theme_set(theme_minimal(base_family = PLOT_FONT))
update_geom_defaults("text", list(family = PLOT_FONT))

# Paths and per-group settings come from the pipeline (--config JSON).
suppressPackageStartupMessages(library(optparse))

.opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character",
              help = "JSON: groups, tumor_samples, seed_files, enumerate_files, peakbed_files, palette, panel_ncol, panel_w, panel_h"),
  make_option("--seed-dir", type = "character", default = "."),
  make_option("--vcf-dir",  type = "character", default = "."),
  make_option("--diff-dir", type = "character", default = "."),
  make_option("--out-dir",  type = "character", default = "."),
  make_option("--genome",   type = "character", help = "Indexed reference FASTA"),
  make_option("--samtools", type = "character", default = "samtools"),
  make_option("--allow-no-on-target", action = "store_true", default = FALSE,
              help = "Allow no on-target sites. Knob for non-differential arm"),
  make_option("--both-strands", action = "store_true", default = FALSE,
              help = "ALT-count plots: also show TS edits below the axis (default: NTS only)")
)))

.cfg <- jsonlite::fromJSON(.opt$config, simplifyVector = FALSE)

seed_dir  <- .opt$`seed-dir`
vcf_dir   <- .opt$`vcf-dir`
diff_dir  <- .opt$`diff-dir`
# spacer sequences are read from the reference (see site_spacer_sequences())
genome_fa <- .opt$genome
SAMTOOLS  <- .opt$samtools
ALLOW_NO_ON_TARGET <- isTRUE(.opt$`allow-no-on-target`)
SHOW_TS <- isTRUE(.opt$`both-strands`)

out_dir <- .opt$`out-dir`
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Stack stage 2's sites alongside stage 1's; `site_source` marks each row's origin.
# The two are not homogeneous -- different mismatch ceilings, different ranking.
INCLUDE_ENUMERATE <- TRUE

# --- per-group config, from --config ---------------------------------------
conditions    <- unlist(.cfg$groups)
# treatment replicates only -- the merged VCF also carries control columns
tumor_samples <- lapply(.cfg$tumor_samples, unlist)
seed_files      <- unlist(.cfg$seed_files)        # stage 1 binding-site TSV per group
enumerate_files <- unlist(.cfg$enumerate_files)   # stage 2 TSV per group
peakbed_files   <- unlist(.cfg$peakbed_files)     # differential peaks BED per group
condition_palette <- unlist(.cfg$palette)
panel_ncol <- .cfg$panel_ncol; panel_w <- .cfg$panel_w; panel_h <- .cfg$panel_h

stopifnot(length(conditions) > 0,
          all(conditions %in% names(tumor_samples)),
          all(conditions %in% names(seed_files)),
          all(conditions %in% names(peakbed_files)),
          all(conditions %in% names(condition_palette)))
cat(sprintf("groups = %s\n", paste(conditions, collapse = ", ")))

snp_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
snp_type_palette <- setNames(
  c("#1565C0", "black", "#C62828", "gray", "#FCAE1E", "#0B6623"),
  snp_types
)

# y half-range: largest stack in either direction, per panel.
panel_ylim <- function(totals_df, cond) {
  lim <- max(abs(totals_df$total[totals_df$condition == cond]), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) 1 else lim
}

# 6-class pyrimidine-collapsed colour class (A>G and T>C merge); strand is
# recovered separately by strand_sign_of().
snp_type_of <- function(ref, alt) {
  key <- paste0(ref, ">", alt)
  collapse <- c(
    "C>A" = "C>A", "G>T" = "C>A",
    "C>G" = "C>G", "G>C" = "C>G",
    "C>T" = "C>T", "G>A" = "C>T",
    "T>A" = "T>A", "A>T" = "T>A",
    "T>C" = "T>C", "A>G" = "T>C",
    "T>G" = "T>G", "A>C" = "T>G"
  )
  unname(collapse[key])
}

COMPLEMENT <- c(A = "T", C = "G", G = "C", T = "A")

# Substitution on the protospacer (NTS) strand: complemented on '-' sites, not
# collapsed (A>G = NTS adenine edit, T>C = TS adenine edit).
ps_substitution <- function(ref, alt, site_strand) {
  neg <- !is.na(site_strand) & site_strand == "-"
  r <- ref
  a <- alt
  r[neg] <- unname(COMPLEMENT[ref[neg]])
  a[neg] <- unname(COMPLEMENT[alt[neg]])
  paste0(r, ">", a)
}

# +1 if the protospacer-strand REF is C/T (matches the collapsed label), else -1.
# The A classes are flipped later (see A_RELABEL).
strand_sign_of <- function(ps_type) {
  ifelse(substr(ps_type, 1, 1) %in% c("C", "T"), 1L, -1L)
}

# Per-replicate FORMAT/AF for one VCF line (AF index looked up per line; NA if uncalled).
sample_afs <- function(format_str, sample_fields) {
  af_idx <- match("AF", strsplit(format_str, ":", fixed = TRUE)[[1]])
  vapply(sample_fields, function(s) {
    # first allele only (multi-allelic rows are dropped later anyway)
    raw <- strsplit(s, ":", fixed = TRUE)[[1]][af_idx]
    if (is.na(raw) || raw == ".") return(NA_real_)
    as.numeric(strsplit(raw, ",", fixed = TRUE)[[1]][1])
  }, numeric(1))
}

# SEM ignoring NAs; NA if fewer than 2 values.
sem <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / sqrt(length(x))
}

# ALT read count (FORMAT/AD alt component) summed across replicates for one VCF line.
sum_alt_ad <- function(format_str, sample_fields) {
  ad_idx <- match("AD", strsplit(format_str, ":", fixed = TRUE)[[1]])
  alt_vals <- vapply(sample_fields, function(s) {
    raw <- strsplit(s, ":", fixed = TRUE)[[1]][ad_idx]
    if (is.na(raw) || raw == ".") return(NA_integer_)
    as.integer(strsplit(raw, ",", fixed = TRUE)[[1]][2])
  }, integer(1))
  sum(alt_vals, na.rm = TRUE)
}

read_variants <- function(vcf_path, samples) {
  lines <- readLines(vcf_path)
  header <- strsplit(lines[startsWith(lines, "#CHROM")], "\t", fixed = TRUE)[[1]]
  sample_col <- match(samples, header)
  if (anyNA(sample_col)) {
    stop(sprintf("sample(s) not found in %s: %s", vcf_path,
                  paste(samples[is.na(sample_col)], collapse = ", ")))
  }

  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) == 0) {
    return(data.frame(chrom = character(0), pos = integer(0), ref = character(0),
                       alt = character(0), af = numeric(0), alt_count = integer(0)))
  }
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  df <- data.frame(
    chrom = vapply(split_lines, `[`, character(1), 1),
    pos   = as.integer(vapply(split_lines, `[`, character(1), 2)),
    ref   = vapply(split_lines, `[`, character(1), 4),
    alt   = vapply(split_lines, `[`, character(1), 5)
  )
  # per-replicate AFs kept for the on-target SEM; `af` is their mean
  df$af_list <- lapply(split_lines, function(fields) {
    unname(sample_afs(fields[9], fields[sample_col]))
  })
  df$af <- vapply(df$af_list, function(v) mean(v, na.rm = TRUE), numeric(1))
  df$alt_count <- vapply(split_lines, function(fields) {
    sum_alt_ad(fields[9], fields[sample_col])
  }, integer(1))
  is_snp <- nchar(df$ref) == 1 & nchar(df$alt) == 1
  n_indel <- sum(!is_snp)
  if (n_indel > 0) {
    cat(sprintf("  %d/%d variants are indels (not a single-base substitution) -- excluded\n",
                n_indel, nrow(df)))
  }
  df <- df[is_snp, ]

  # Mutect2 can emit the same SNV twice (two haplotype events); keep the
  # highest-AF record so an edit isn't counted twice. Row order is preserved.
  key <- paste(df$chrom, df$pos, df$ref, df$alt, sep = ":")
  dup_keys <- unique(key[duplicated(key)])
  if (length(dup_keys) > 0) {
    drop <- logical(nrow(df))
    for (k in dup_keys) {
      idx <- which(key == k)
      af <- df$af[idx]
      af[is.na(af)] <- -Inf              # an uncalled AF never wins
      drop[setdiff(idx, idx[which.max(af)])] <- TRUE
    }
    cat(sprintf("  %d duplicate SNV record(s) across %d position(s) collapsed (highest AF kept)\n",
                sum(drop), length(dup_keys)))
    df <- df[!drop, ]
  }

  df$snp_type <- snp_type_of(df$ref, df$alt)
  df
}

# Step 7 (enumerate) sites with the 3 nt PAM trimmed (strand-aware) to match step 6's
# spacer-only intervals. NULL with a warning if the file is missing.
PAM_LEN <- 3L

read_enumerate_sites <- function(cond) {
  p <- file.path(seed_dir, enumerate_files[[cond]])
  if (!file.exists(p)) {
    warning(sprintf("%s: no %s -- step 7 not run for this condition, using step 6 sites only",
                     cond, basename(p)))
    return(NULL)
  }
  df <- read.table(p, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df <- df[df$start != ".", ]
  if (nrow(df) == 0) return(NULL)
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  plus <- df$strand == "+"
  df$end[plus]    <- df$end[plus]    - PAM_LEN   # PAM at the high end
  df$start[!plus] <- df$start[!plus] + PAM_LEN   # PAM at the low end
  df$seed_hits <- 0L
  df
}

# Seed + enumerate sites for one condition, stacked; `site_source` records the origin.
read_sites <- function(tsv_path, cond = NULL) {
  df <- read.table(tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df <- df[df$start != ".", ]
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  df$site_source <- "seed"
  if (INCLUDE_ENUMERATE && !is.null(cond)) {
    enum <- read_enumerate_sites(cond)
    if (!is.null(enum)) {
      enum$site_source <- "enumerate"
      # step 7 only covers peaks step 6 left unscored, so no double-counting
      df <- rbind(df[, intersect(names(df), names(enum))],
                   enum[, intersect(names(df), names(enum))])
    }
  }
  # on-target = the unique 0-mismatch/0-bulge hit (needs SUMMIT_WIDTH=0 in guidescan)
  df$on_target <- as.integer(df$mismatches) == 0 &
    as.integer(df$rna_bulges) == 0 & as.integer(df$dna_bulges) == 0
  df[, c("chrom", "start", "end", "strand", "on_target", "site_source",
         "peak_chrom", "peak_start", "peak_end")]
}

position_in_protospacer <- function(vcf_path, tsv_path, samples, cond = NULL) {
  variants <- read_variants(vcf_path, samples)
  sites <- read_sites(tsv_path, cond)

  n_on_target <- sum(sites$on_target)
  if (n_on_target > 1 || (n_on_target == 0 && !ALLOW_NO_ON_TARGET)) {
    stop(sprintf("%s: expected exactly 1 on-target (0 mismatch, 0 bulge) site, found %d",
                  tsv_path, n_on_target))
  }

  # sites are 0-based half-open, spacer only (no PAM) -> 1-based closed, positions 1..20
  sites_gr <- GRanges(seqnames = sites$chrom,
                       ranges = IRanges(start = sites$start + 1L, end = sites$end),
                       strand = sites$strand)
  variants_gr <- GRanges(seqnames = variants$chrom,
                          ranges = IRanges(start = variants$pos, width = 1L))

  hit <- findOverlaps(variants_gr, sites_gr, select = "first")
  n_total <- length(variants_gr)
  n_matched <- sum(!is.na(hit))
  cat(sprintf("  %d/%d substitutions (%s) matched a scored protospacer\n",
              n_matched, n_total,
              if (n_total > 0) sprintf("%.1f%%", 100 * n_matched / n_total) else "n/a"))

  matched <- !is.na(hit)
  site_start <- start(sites_gr)[hit[matched]]
  site_end <- end(sites_gr)[hit[matched]]
  site_strand <- as.character(strand(sites_gr)[hit[matched]])
  variant_pos <- start(variants_gr)[matched]

  position <- ifelse(site_strand == "+",
                      variant_pos - site_start + 1L,
                      site_end - variant_pos + 1L)
  # ps_type needs the matched site's strand, so it's resolved here rather than in read_variants()
  ps_type <- ps_substitution(variants$ref[matched], variants$alt[matched], site_strand)

  result <- data.frame(position = position, snp_type = variants$snp_type[matched],
                        ps_type = ps_type,
                        strand_sign = strand_sign_of(ps_type),
                        site_strand = site_strand,
                        site_source = sites$site_source[hit[matched]],
                        site_id = hit[matched], af = variants$af[matched],
                        alt_count = variants$alt_count[matched],
                        chrom = variants$chrom[matched], genomic_pos = variants$pos[matched],
                        ref = variants$ref[matched], alt = variants$alt[matched],
                        on_target = sites$on_target[hit[matched]],
                        peak_id = paste(sites$peak_chrom, sites$peak_start,
                                        sites$peak_end)[hit[matched]])
  result$af_list <- variants$af_list[matched]
  result
}

# Every matched substitution per condition, on-target included (plots filter it out).
hits_by_condition <- list()
for (cond in conditions) {
  cat(sprintf("%s:\n", cond))
  vcf_path <- file.path(vcf_dir, cond, sprintf("%s_vs_control_relcontrol_filt_no-alt-chr.vcf", cond))
  tsv_path <- file.path(seed_dir, seed_files[[cond]])
  hits_by_condition[[cond]] <- position_in_protospacer(vcf_path, tsv_path, tumor_samples[[cond]], cond)
}

# On-target report: mean AF +/- SEM across replicates for every SNP in the
# on-target spacer (may be more than one per condition).
on_target_summary <- do.call(rbind, lapply(conditions, function(cond) {
  ot <- hits_by_condition[[cond]]
  ot <- ot[ot$on_target, ]
  if (nrow(ot) == 0) {
    warning(sprintf("%s: no on-target SNP found in the on-target window", cond))
    return(data.frame(condition = cond, chrom = NA_character_, genomic_pos = NA_integer_,
                       ref = NA_character_, alt = NA_character_,
                       n_replicates = NA_integer_, mean_af = NA_real_, sem_af = NA_real_,
                       af_by_replicate = NA_character_))
  }
  do.call(rbind, lapply(seq_len(nrow(ot)), function(i) {
    vals <- ot$af_list[[i]]
    valid <- vals[!is.na(vals)]
    # af_by_replicate shows which replicate carried each call (vals is in tumor_samples order)
    data.frame(condition = cond, chrom = ot$chrom[i], genomic_pos = ot$genomic_pos[i],
               ref = ot$ref[i], alt = ot$alt[i], n_replicates = length(valid),
               mean_af = if (length(valid) > 0) mean(valid) else NA_real_,
               sem_af = sem(vals),
               af_by_replicate = paste(
                 sprintf("%s=%s", tumor_samples[[cond]],
                          vapply(vals, function(v) if (is.na(v)) "NA" else as.character(v),
                                  character(1))),
                 collapse = ";"))
  }))
}))
write.csv(on_target_summary, file.path(out_dir, "on_target_af.csv"), row.names = FALSE)

# Panel title per condition, using the highest-AF on-target SNP.
on_target_label <- setNames(vapply(conditions, function(cond) {
  ot <- on_target_summary[on_target_summary$condition == cond, ]
  ot <- ot[!is.na(ot$mean_af), ]
  if (nrow(ot) == 0) return(sprintf("%s (on-target: not found)", cond))
  best <- ot[which.max(ot$mean_af), ]
  sprintf("%s (on-target max AF %.2f\u00b1%.2f)", cond, best$mean_af, best$sem_af)
}, character(1)), conditions)

# add subtitle
RLOOP_SIGNATURES <- c("A>G", "C>T")

rloop_label <- setNames(vapply(conditions, function(cond) {
  hits <- hits_by_condition[[cond]]
  n_peaks <- length(readLines(file.path(diff_dir, peakbed_files[[cond]]))) -
    length(unique(hits$peak_id[hits$on_target]))
  hits <- hits[!hits$on_target & !is.na(hits$af), ]
  paste(vapply(RLOOP_SIGNATURES, function(sig) {
    is_sig <- hits$ps_type == sig
    if (!any(is_sig) || n_peaks <= 0) return(sprintf("%s: no edited R-loops", sig))
    per_peak <- tapply(hits$af[is_sig], hits$peak_id[is_sig], max)
    sprintf("%s: %.1f%% of R-loops edited (%s/%s), median AF %.1f%%",
            sig, 100 * length(per_peak) / n_peaks,
            formatC(length(per_peak), format = "d", big.mark = ","),
            formatC(n_peaks, format = "d", big.mark = ","),
            100 * median(per_peak))
  }, character(1)), collapse = "\n")
}, character(1)), conditions)

write.csv(data.frame(condition = conditions,
                     rloop_summary = gsub("\n", " | ", rloop_label[conditions])),
          file.path(out_dir, "rloop_offtarget_summary.csv"), row.names = FALSE)

# Stack height per position, summed separately for up and down so they don't cancel.
signed_totals <- function(df, value_col) {
  df %>%
    mutate(dir = ifelse(.data[[value_col]] < 0, "down", "up")) %>%
    group_by(condition, position, dir) %>%
    summarise(total = sum(.data[[value_col]]), .groups = "drop")
}

# Legend label per class with its NTS (up) / TS (down) totals; NTS only unless show_ts.
signed_type_labels <- function(df, value_col, fmt = "%s", show_ts = TRUE) {
  v <- df[[value_col]]
  agg <- df %>%
    mutate(.v = v) %>%
    group_by(snp_type) %>%
    summarise(up = sum(.v[.v > 0]), down = -sum(.v[.v < 0]), .groups = "drop")
  if (!show_ts) {
    return(setNames(sprintf(paste0("%s (NTS ", fmt, ")"), agg$snp_type, agg$up), agg$snp_type))
  }
  setNames(sprintf(paste0("%s (NTS ", fmt, " / TS ", fmt, ")"),
                    agg$snp_type, agg$up, agg$down), agg$snp_type)
}

x_range <- range(unlist(lapply(hits_by_condition, function(h) h$position[!h$on_target])))

# kept short: repeated in every panel
X_LAB <- "Position in protospacer (1 = PAM-distal)"

# ALT reads (summed over replicates and off-target sites) per position x type x
# strand, zero-filled. Up = NTS, down = TS.
all_altcount <- do.call(rbind, lapply(conditions, function(cond) {
  hits <- hits_by_condition[[cond]]
  hits <- hits[!hits$on_target, ]
  hits$snp_type <- factor(hits$snp_type, levels = snp_types)
  if (nrow(hits) == 0) return(NULL)
  positions <- sort(unique(hits$position))
  grid <- expand.grid(position = positions,
                       snp_type = factor(snp_types, levels = snp_types),
                       strand_sign = c(-1L, 1L))
  altcount_by_group <- hits %>% group_by(position, snp_type, strand_sign) %>%
    summarise(alt_count = sum(alt_count), .groups = "drop")
  grid %>%
    left_join(altcount_by_group, by = c("position", "snp_type", "strand_sign")) %>%
    mutate(alt_count = ifelse(is.na(alt_count), 0, alt_count) * strand_sign, condition = cond)
}))
all_altcount$condition <- factor(all_altcount$condition, levels = conditions)

# Relabel T>A/T>C/T>G as A>T/A>G/A>C and flip their sign so up = A>X on the NTS.
# Done before write.csv so the CSV matches the figure; the CSV keeps both strands
# even when the figure is NTS only.
A_RELABEL <- c("T>A" = "A>T", "T>C" = "A>G", "T>G" = "A>C")
plot_snp_types <- c("C>A", "C>G", "C>T", "A>T", "A>G", "A>C")
plot_snp_palette <- setNames(unname(snp_type_palette[snp_types]), plot_snp_types)

flip <- as.character(all_altcount$snp_type) %in% names(A_RELABEL)
all_altcount$alt_count[flip]   <- -all_altcount$alt_count[flip]
all_altcount$strand_sign[flip] <- -all_altcount$strand_sign[flip]
relabeled <- as.character(all_altcount$snp_type)
relabeled[flip] <- unname(A_RELABEL[relabeled[flip]])
all_altcount$snp_type <- factor(relabeled, levels = plot_snp_types)

write.csv(all_altcount, file.path(out_dir, "variant_altcount_in_protospacer.csv"), row.names = FALSE)

# NTS (strand_sign = +1) only unless --both-strands.
plot_altcount <- if (SHOW_TS) all_altcount else all_altcount[all_altcount$strand_sign == 1L, ]
altcount_position_totals <- signed_totals(plot_altcount, "alt_count")

# One ggplot per condition (not facet_wrap) so each legend shows its own totals.
altcount_plots <- lapply(conditions, function(cond) {
  df <- plot_altcount[plot_altcount$condition == cond, ]
  type_labels <- signed_type_labels(df, "alt_count", "%d", show_ts = SHOW_TS)
  altcount_y_lim <- panel_ylim(altcount_position_totals, cond)

  ggplot(df, aes(x = position, y = alt_count, fill = snp_type)) +
    geom_col(color = "black", linewidth = 0.2, width = 0.8, alpha = 0.9) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
    # strand labels pinned to the panel's top-left / bottom-left corners
    annotate("text", x = -Inf, y = Inf, label = "NTS", hjust = -0.3, vjust = 1.5,
             fontface = "bold", size = 3.5, family = PLOT_FONT) +
    (if (SHOW_TS) annotate("text", x = -Inf, y = -Inf, label = "TS", hjust = -0.3, vjust = -0.5,
                           fontface = "bold", size = 3.5, family = PLOT_FONT)) +
    scale_fill_manual(values = plot_snp_palette, labels = type_labels, name = "SNP type", drop = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = function(v) format(abs(v), big.mark = ",")) +
    coord_cartesian(xlim = x_range,
                    ylim = c(if (SHOW_TS) -altcount_y_lim else 0, altcount_y_lim)) +
    labs(x = X_LAB, y = "Total ALT reads (3 reps)", title = on_target_label[[cond]],
         subtitle = rloop_label[[cond]]) +
    theme_minimal(base_size = 11, base_family = PLOT_FONT) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", color = condition_palette[[cond]]),
          plot.subtitle = element_text(size = 7, lineheight = 1.2, color = "grey25"))
})

p_altcount <- wrap_plots(altcount_plots, ncol = panel_ncol) +
  plot_annotation(theme = theme(text = element_text(family = PLOT_FONT)))

ggsave(file.path(out_dir, "variant_altcount_in_protospacer.pdf"), p_altcount, width = panel_w, height = panel_h)
ggsave(file.path(out_dir, "variant_altcount_in_protospacer.png"), p_altcount, width = panel_w, height = panel_h, dpi = 200)

# Per-edit weight for AF-normalized edit counts: mean AF over all treatment replicates,
# an uncalled replicate counting as AF 0. A weighted sum equals the edit count if every
# edit were at AF 1 in every replicate.
mean_af_all_reps <- function(af_list) {
  vapply(af_list, function(v) sum(v, na.rm = TRUE) / length(v), numeric(1))
}

# Differential peaks (BED line count) per AF-weighted off-target NTS edit, one set of
# outputs per signature: peaks_per_edit_AG.* (ABE) and peaks_per_edit_CT.* (CBE).
peak_counts <- data.frame(
  condition = conditions,
  n_peaks = vapply(conditions, function(cond) {
    length(readLines(file.path(diff_dir, peakbed_files[[cond]])))
  }, integer(1))
)
peak_counts$condition <- factor(peak_counts$condition, levels = conditions)

# `collapsed` = the 6-class label of the signature; its count (both strands) is
# kept as n_collapsed_TC / n_collapsed_CT for comparison with older runs.
PEAKS_PER_EDIT_SIGNATURES <- list(
  AG = list(sig = "A>G", collapsed = "T>C"),
  CT = list(sig = "C>T", collapsed = "C>T")
)

for (tag in names(PEAKS_PER_EDIT_SIGNATURES)) {
  s <- PEAKS_PER_EDIT_SIGNATURES[[tag]]
  collapsed_col <- paste0("n_collapsed_", sub(">", "", s$collapsed, fixed = TRUE))

  signature_edit_counts <- do.call(rbind, lapply(conditions, function(cond) {
    hits <- hits_by_condition[[cond]]
    hits <- hits[!hits$on_target, ]
    is_sig <- hits$ps_type == s$sig
    df <- data.frame(condition = cond,
                     signature = sprintf("%s (protospacer strand)", s$sig),
                     n_signature_edits = sum(is_sig),
                     n_signature_edits_af = sum(mean_af_all_reps(hits$af_list[is_sig])))
    df[[collapsed_col]] <- sum(hits$snp_type == s$collapsed)
    df
  }))

  peaks_per_edit <- merge(peak_counts, signature_edit_counts, by = "condition")
  peaks_per_edit$condition <- factor(peaks_per_edit$condition, levels = conditions)
  peaks_per_edit <- peaks_per_edit[order(peaks_per_edit$condition), ]
  # denominator is AF-weighted; peaks_per_edit_count keeps the unweighted ratio
  peaks_per_edit$peaks_per_edit <- ifelse(peaks_per_edit$n_signature_edits_af > 0,
                                           peaks_per_edit$n_peaks / peaks_per_edit$n_signature_edits_af, NA)
  peaks_per_edit$peaks_per_edit_count <- ifelse(peaks_per_edit$n_signature_edits > 0,
                                                 peaks_per_edit$n_peaks / peaks_per_edit$n_signature_edits, NA)

  write.csv(peaks_per_edit, file.path(out_dir, sprintf("peaks_per_edit_%s.csv", tag)), row.names = FALSE)

  p_ratio <- ggplot(peaks_per_edit, aes(x = condition, y = peaks_per_edit, fill = condition)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.7, alpha = 0.9) +
    # ratio / peaks / AF-weighted edits, one per line to stay within the bar
    geom_text(aes(label = sprintf("%s\n%s peaks\n%s AF-wt. edits",
                                    ifelse(is.na(peaks_per_edit), "NA",
                                           formatC(peaks_per_edit, format = "f", digits = 1, big.mark = ",")),
                                    formatC(n_peaks, format = "d", big.mark = ","),
                                    formatC(n_signature_edits_af, format = "f", digits = 1, big.mark = ","))),
              vjust = -0.2, size = 2.7, lineheight = 0.95) +
    scale_fill_manual(values = condition_palette) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.3))) +
    labs(x = NULL,
         y = sprintf("Differential peaks per AF-weighted %s edit (NTS)", s$sig)) +
    theme_minimal(base_size = 12, base_family = PLOT_FONT) +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank())

  ggsave(file.path(out_dir, sprintf("peaks_per_edit_%s.pdf", tag)), p_ratio, width = 8.5, height = 6)
  ggsave(file.path(out_dir, sprintf("peaks_per_edit_%s.png", tag)), p_ratio, width = 8.5, height = 6, dpi = 200)
}

# peak_snp_summary.csv: one row per differential peak (on-target included), with
# mean AF (af_1..20) and NTS-oriented substitution (type_1..20) at each position.
PROTOSPACER_LEN <- 20L

# Protospacer-strand sequence of each site from the reference (one batched faidx);
# '-' sites are reverse-complemented so index 1 is PAM-distal.
site_spacer_sequences <- function(sites) {
  if (nrow(sites) == 0) return(character(0))
  regfile <- tempfile(fileext = ".txt")
  # site_start is 0-based half-open; faidx wants 1-based inclusive
  writeLines(sprintf("%s:%d-%d", sites$site_chrom, sites$site_start + 1L, sites$site_end),
             regfile)
  raw <- system2(SAMTOOLS, c("faidx", shQuote(genome_fa), "-r", shQuote(regfile)),
                  stdout = TRUE, stderr = FALSE)
  unlink(regfile)
  # collapse the multi-FASTA back to one string per requested region, in order
  idx <- cumsum(startsWith(raw, ">"))
  seqs <- vapply(split(raw[!startsWith(raw, ">")], idx[!startsWith(raw, ">")]),
                  function(x) toupper(paste0(x, collapse = "")), character(1))
  if (length(seqs) != nrow(sites)) {
    stop(sprintf("samtools returned %d sequence(s) for %d site(s)", length(seqs), nrow(sites)))
  }
  neg <- !is.na(sites$strand) & sites$strand == "-"
  seqs[neg] <- vapply(seqs[neg], function(x)
    paste(rev(strsplit(chartr("ACGT", "TGCA", x), "", fixed = TRUE)[[1]]), collapse = ""),
    character(1))
  unname(seqs)
}

# Every row of the guidescan TSV (one per peak), unscored ones included.
read_sites_all <- function(tsv_path, cond = NULL) {
  df <- read.table(tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = "character")
  num <- function(x) suppressWarnings(as.numeric(ifelse(x == ".", NA, x)))
  int <- function(x) suppressWarnings(as.integer(ifelse(x == ".", NA, x)))
  out <- data.frame(
    peak_chrom      = df$peak_chrom,
    peak_start      = int(df$peak_start),
    peak_end        = int(df$peak_end),
    guidescan_match = df$start != ".",
    seed_hits       = int(df$seed_hits),
    site_chrom      = ifelse(df$start == ".", NA_character_, df$chrom),
    site_start      = int(df$start),
    site_end        = int(df$end),
    strand          = ifelse(df$strand == ".", NA_character_, df$strand),
    mismatches      = int(df$mismatches),
    rna_bulges      = int(df$rna_bulges),
    dna_bulges      = int(df$dna_bulges),
    cfd             = num(df$cfd),
    # protospacer orientation, 20 nt spacer + 3 nt PAM
    match_sequence  = ifelse(df$start == ".", NA_character_, df$match_sequence),
    site_source     = ifelse(df$start == ".", NA_character_, "seed"),
    stringsAsFactors = FALSE
  )

  # fill step 6's unscored rows from step 7 (filled, not appended, to avoid duplicate peaks)
  if (INCLUDE_ENUMERATE && !is.null(cond)) {
    enum <- read_enumerate_sites(cond)
    if (!is.null(enum)) {
      key  <- paste(out$peak_chrom, out$peak_start, out$peak_end)
      ekey <- paste(enum$peak_chrom, enum$peak_start, enum$peak_end)
      i  <- match(ekey, key)
      ok <- !is.na(i) & !out$guidescan_match[i]   # only ever fill a blank row
      i  <- i[ok]
      e  <- enum[ok, ]
      out$guidescan_match[i] <- TRUE
      out$site_chrom[i]      <- e$chrom
      out$site_start[i]      <- e$start
      out$site_end[i]        <- e$end
      out$strand[i]          <- e$strand
      out$mismatches[i]      <- as.integer(e$mismatches)
      out$rna_bulges[i]      <- as.integer(e$rna_bulges)
      out$dna_bulges[i]      <- as.integer(e$dna_bulges)
      out$cfd[i]             <- as.numeric(e$cfd)
      out$match_sequence[i]  <- e$match_sequence
      out$site_source[i]     <- "enumerate"
    }
  }
  out
}

peak_snp_table <- function(cond, vcf_path, tsv_path, samples) {
  peaks <- read_sites_all(tsv_path, cond)
  peaks$on_target <- !is.na(peaks$mismatches) & peaks$mismatches == 0L &
    peaks$rna_bulges == 0L & peaks$dna_bulges == 0L
  peaks$n_snps <- 0L

  af_mat <- matrix(NA_real_, nrow(peaks), PROTOSPACER_LEN)
  type_mat <- matrix(NA_character_, nrow(peaks), PROTOSPACER_LEN)

  matched_idx <- which(peaks$guidescan_match)
  if (length(matched_idx) > 0) {
    variants <- read_variants(vcf_path, samples)
    if (nrow(variants) > 0) {
      sites_gr <- GRanges(seqnames = peaks$site_chrom[matched_idx],
                          ranges = IRanges(start = peaks$site_start[matched_idx] + 1L,
                                            end = peaks$site_end[matched_idx]),
                          strand = peaks$strand[matched_idx])
      variants_gr <- GRanges(seqnames = variants$chrom,
                             ranges = IRanges(start = variants$pos, width = 1L))
      hit <- findOverlaps(variants_gr, sites_gr, select = "first")
      keep <- which(!is.na(hit))
      for (v in keep) {
        row <- matched_idx[hit[v]]          # back to the full-peak table
        pos <- if (as.character(strand(sites_gr)[hit[v]]) == "+") {
          variants$pos[v] - start(sites_gr)[hit[v]] + 1L
        } else {
          end(sites_gr)[hit[v]] - variants$pos[v] + 1L
        }
        peaks$n_snps[row] <- peaks$n_snps[row] + 1L
        if (pos < 1L || pos > PROTOSPACER_LEN) next   # defensive; shouldn't occur
        # if two SNPs share a position, keep the higher-AF one
        if (is.na(af_mat[row, pos]) ||
            (!is.na(variants$af[v]) && variants$af[v] > af_mat[row, pos])) {
          af_mat[row, pos] <- variants$af[v]
          type_mat[row, pos] <- ps_substitution(variants$ref[v], variants$alt[v], peaks$strand[row])
        }
      }
    }
  }

  colnames(af_mat) <- sprintf("af_%d", seq_len(PROTOSPACER_LEN))
  colnames(type_mat) <- sprintf("type_%d", seq_len(PROTOSPACER_LEN))
  cbind(data.frame(condition = cond, stringsAsFactors = FALSE), peaks,
        as.data.frame(af_mat), as.data.frame(type_mat, stringsAsFactors = FALSE))
}

peak_snp_summary <- do.call(rbind, lapply(conditions, function(cond) {
  vcf_path <- file.path(vcf_dir, cond, sprintf("%s_vs_control_relcontrol_filt_no-alt-chr.vcf", cond))
  tsv_path <- file.path(seed_dir, seed_files[[cond]])
  tab <- peak_snp_table(cond, vcf_path, tsv_path, tumor_samples[[cond]])
  cat(sprintf("%s: %d peaks, %d with a guidescan match, %d with >=1 SNP in the protospacer\n",
              cond, nrow(tab), sum(tab$guidescan_match), sum(tab$n_snps > 0)))
  tab
}))

write.csv(peak_snp_summary, file.path(out_dir, "peak_snp_summary.csv"), row.names = FALSE)
cat(sprintf("Wrote %s: %d rows x %d columns\n",
            file.path(out_dir, "peak_snp_summary.csv"),
            nrow(peak_snp_summary), ncol(peak_snp_summary)))

# Substrate-normalized rate: off-target NTS signature edits per available NTS substrate
# base, one set of outputs per signature: substrate_normalized_summary_AG.* (A>G per
# adenine) and substrate_normalized_summary_CT.* (C>T per cytosine).
# Clean within a locus; across loci, seed/guide composition affects the denominator.
substrate_by_position <- do.call(rbind, lapply(conditions, function(cond) {
  sites <- read_sites_all(file.path(seed_dir, seed_files[[cond]]), cond)
  sites <- sites[sites$guidescan_match, ]
  on_t <- !is.na(sites$mismatches) & sites$mismatches == 0L &
    sites$rna_bulges == 0L & sites$dna_bulges == 0L
  sites <- sites[!on_t, ]
  spacer <- site_spacer_sequences(sites)
  spacer <- spacer[nchar(spacer) > 0]
  # bulged sites count only at the positions they span (NA tail / truncated at 20)
  chars <- matrix(NA_character_, length(spacer), PROTOSPACER_LEN)
  for (i in seq_along(spacer)) {
    cs <- strsplit(spacer[i], "", fixed = TRUE)[[1]]
    k <- min(length(cs), PROTOSPACER_LEN)
    chars[i, seq_len(k)] <- cs[seq_len(k)]
  }
  data.frame(condition = cond, position = seq_len(PROTOSPACER_LEN),
             n_sites = length(spacer),
             n_A = as.integer(colSums(chars == "A", na.rm = TRUE)),
             n_C = as.integer(colSums(chars == "C", na.rm = TRUE)),
             n_G = as.integer(colSums(chars == "G", na.rm = TRUE)),
             n_T = as.integer(colSums(chars == "T", na.rm = TRUE)),
             stringsAsFactors = FALSE)
}))

# observed strand-resolved edits per position, off-target only. n_TC / n_GA are the
# same substitutions on the TS (seen from the NTS), kept for comparison.
# *_af columns weight each edit by its AF instead of 1:
#   _af        mean_af_all_reps() (uncalled replicate = AF 0), as in peaks_per_edit
#   _af_called mean AF over the replicates that called it (the `af` used elsewhere)
edits_by_position <- do.call(rbind, lapply(conditions, function(cond) {
  hits <- hits_by_condition[[cond]]
  hits <- hits[!hits$on_target, ]
  af_all <- mean_af_all_reps(hits$af_list)
  af_called <- ifelse(is.na(hits$af), 0, hits$af)
  count <- function(type) as.integer(tabulate(hits$position[hits$ps_type == type], PROTOSPACER_LEN))
  weighted <- function(type, w) {
    sel <- hits$ps_type == type
    vapply(seq_len(PROTOSPACER_LEN), function(p) sum(w[sel & hits$position == p]), numeric(1))
  }
  data.frame(condition = cond, position = seq_len(PROTOSPACER_LEN),
             n_AG = count("A>G"), n_TC = count("T>C"),
             n_CT = count("C>T"), n_GA = count("G>A"),
             n_AG_af = weighted("A>G", af_all), n_AG_af_called = weighted("A>G", af_called),
             n_CT_af = weighted("C>T", af_all), n_CT_af_called = weighted("C>T", af_called),
             stringsAsFactors = FALSE)
}))

substrate_rate <- merge(substrate_by_position, edits_by_position,
                         by = c("condition", "position"))
substrate_rate$condition <- factor(substrate_rate$condition, levels = conditions)

# Denominator: protospacer-strand substrate bases (`sub`). `comp` / comp edit = the
# complementary base and substitution, reported alongside but not used in the rate.
SUBSTRATE_SIGNATURES <- list(
  AG = list(sub = "A", edit = "G", substrate = "adenines"),
  CT = list(sub = "C", edit = "T", substrate = "cytosines")
)

# ---- one normalized figure per signature, vs. the un-normalized ratio ------
for (tag in names(SUBSTRATE_SIGNATURES)) {
  s <- SUBSTRATE_SIGNATURES[[tag]]
  comp <- COMPLEMENT[[s$sub]]
  comp_edit <- COMPLEMENT[[s$edit]]
  sig_col  <- paste0("n_", s$sub, s$edit)
  comp_col <- paste0("n_", comp, comp_edit)

  substrate_summary <- substrate_rate %>%
    group_by(condition) %>%
    summarise(n_sites = first(n_sites),
              available_sub = sum(.data[[paste0("n_", s$sub)]]),
              available_comp = sum(.data[[paste0("n_", comp)]]),
              edits_sig = sum(.data[[sig_col]]),
              edits_sig_af = sum(.data[[paste0(sig_col, "_af")]]),
              edits_sig_af_called = sum(.data[[paste0(sig_col, "_af_called")]]),
              edits_comp = sum(.data[[comp_col]]),
              .groups = "drop") %>%
    # rates use the AF-weighted numerator (equals the edit count if every edit had AF 1)
    mutate(available_denom = available_sub,
           mean_denom_per_site = available_denom / n_sites,
           edits_per_10k_substrate = 10000 * edits_sig_af / available_denom,
           edits_per_10k_substrate_count = 10000 * edits_sig / available_denom,
           # the un-normalized metric, for side-by-side comparison
           edits_per_site = edits_sig_af / n_sites) %>%
    select(condition, n_sites, available_sub, available_comp, available_denom,
           edits_sig, edits_sig_af, edits_sig_af_called, edits_comp, mean_denom_per_site,
           edits_per_10k_substrate, edits_per_10k_substrate_count, edits_per_site)
  # e.g. available_A / available_T / edits_AG(_af, _af_called) / edits_TC for A>G
  sig_name <- paste0("edits_", s$sub, s$edit)
  names(substrate_summary) <- sub("^available_sub$", paste0("available_", s$sub),
                              sub("^available_comp$", paste0("available_", comp),
                              sub("^edits_sig", sig_name,
                              sub("^edits_comp$", paste0("edits_", comp, comp_edit),
                                  names(substrate_summary)))))
  substrate_summary$condition <- factor(substrate_summary$condition, levels = conditions)
  edits_col <- paste0(sig_name, "_af")

  write.csv(substrate_summary, file.path(out_dir, sprintf("substrate_normalized_summary_%s.csv", tag)),
            row.names = FALSE)

  plot_df <- substrate_summary
  plot_df$edits <- plot_df[[edits_col]]
  p_norm <- ggplot(plot_df, aes(x = condition, y = edits_per_10k_substrate, fill = condition)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.7, alpha = 0.9) +
    # rate / available bases / AF-weighted edits, one per line to stay within the bar
    geom_text(aes(label = sprintf("%.2f\n%s %s\n%s AF-wt. edits", edits_per_10k_substrate,
                                    formatC(available_denom, format = "d", big.mark = ","), s$sub,
                                    formatC(edits, format = "f", digits = 1, big.mark = ","))),
              vjust = -0.2, size = 2.9, lineheight = 0.95) +
    scale_fill_manual(values = condition_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(x = NULL,
         y = sprintf("AF-weighted %s>%s edits per 10,000 available %s", s$sub, s$edit, s$substrate)) +
    theme_minimal(base_size = 12, base_family = PLOT_FONT) +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank())

  ggsave(file.path(out_dir, sprintf("substrate_normalized_summary_%s.pdf", tag)), p_norm, width = 8, height = 6)
  ggsave(file.path(out_dir, sprintf("substrate_normalized_summary_%s.png", tag)), p_norm, width = 8, height = 6, dpi = 200)

  cat(sprintf("\nSubstrate-normalized summary, %s>%s (off-target only):\n", s$sub, s$edit))
  for (i in seq_len(nrow(substrate_summary))) {
    r <- substrate_summary[i, ]
    cat(sprintf("  %-14s sites=%5d  substrate=%6d (%.2f/site; %s=%d %s=%d)  edits=%4d (AF-weighted %.1f)  %7.2f per 10k %s   [unnormalized: %.4f edits/site]\n",
                as.character(r$condition), r$n_sites, r$available_denom, r$mean_denom_per_site,
                s$sub, r[[paste0("available_", s$sub)]], comp, r[[paste0("available_", comp)]],
                r[[sig_name]], r[[edits_col]], r$edits_per_10k_substrate, s$sub, r$edits_per_site))
  }
}


# VCFs per condition x editor: observed signature edits (numerator) and every possible
# substrate base (denominator). REF/ALT are reference-strand; INFO/PS_* give the NTS view.
vcf_out_dir <- file.path(out_dir, "vcf")
dir.create(vcf_out_dir, showWarnings = FALSE, recursive = TRUE)

# reference contig order, for sorting and for the possible file's header
fai <- read.table(paste0(genome_fa, ".fai"), sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("contig", "length", "offset", "linebases", "linewidth"))
editors <- c("ABE", "CBE")

EDITOR_SPECS <- list(
  ABE = list(tag = "AG", sub = "A", edit = "G", substrate = "adenine",
             deamination = "adenine deamination",
             name = "ABE (adenine base editor; ABE8e in this experiment)"),
  CBE = list(tag = "CT", sub = "C", edit = "T", substrate = "cytosine",
             deamination = "cytosine deamination",
             name = "CBE (cytosine base editor)")
)

read_variants_raw <- function(vcf_path, samples) {
  lines <- readLines(vcf_path)
  meta <- lines[startsWith(lines, "##")]
  chrom_line <- lines[startsWith(lines, "#CHROM")]
  header <- strsplit(chrom_line, "\t", fixed = TRUE)[[1]]
  sample_col <- match(samples, header)
  if (anyNA(sample_col)) {
    stop(sprintf("sample(s) not found in %s: %s", vcf_path,
                 paste(samples[is.na(sample_col)], collapse = ", ")))
  }

  body <- lines[!startsWith(lines, "#")]
  empty <- data.frame(chrom = character(0), pos = integer(0), id = character(0),
                      ref = character(0), alt = character(0), qual = character(0),
                      filter = character(0), fmt = character(0),
                      af = numeric(0), stringsAsFactors = FALSE)
  if (length(body) == 0) return(list(df = empty, meta = meta))

  f <- strsplit(body, "\t", fixed = TRUE)
  g <- function(i) vapply(f, `[`, character(1), i)
  df <- data.frame(chrom = g(1), pos = as.integer(g(2)), id = g(3),
                   ref = g(4), alt = g(5), qual = g(6), filter = g(7), fmt = g(9),
                   stringsAsFactors = FALSE)
  df$smp <- lapply(f, function(x) x[sample_col])
  df$af  <- vapply(f, function(x) mean(unname(sample_afs(x[9], x[sample_col])), na.rm = TRUE),
                   numeric(1))

  is_snp <- nchar(df$ref) == 1 & nchar(df$alt) == 1
  if (sum(!is_snp) > 0) {
    cat(sprintf("  %d/%d variants are indels -- excluded\n", sum(!is_snp), nrow(df)))
  }
  df <- df[is_snp, ]

  key <- paste(df$chrom, df$pos, df$ref, df$alt, sep = ":")
  dup_keys <- unique(key[duplicated(key)])
  if (length(dup_keys) > 0) {
    drop <- logical(nrow(df))
    for (k in dup_keys) {
      idx <- which(key == k)
      af <- df$af[idx]
      af[is.na(af)] <- -Inf          # an uncalled AF never wins
      drop[setdiff(idx, idx[which.max(af)])] <- TRUE
    }
    cat(sprintf("  %d duplicate SNV record(s) across %d position(s) collapsed (highest AF kept)\n",
                sum(drop), length(dup_keys)))
    df <- df[!drop, ]
  }
  list(df = df, meta = meta)
}

info_lines_common <- function(spec) c(
  '##INFO=<ID=PS_COND,Number=1,Type=String,Description="Experimental condition (differential peak set) this record came from">',
  sprintf('##INFO=<ID=PS_EDITOR,Number=1,Type=String,Description="Base editor this file describes: %s. Its protospacer-strand signature is %s>%s">',
          spec$name, spec$sub, spec$edit),
  '##INFO=<ID=PS_POS,Number=1,Type=Integer,Description="1-based position within the 20 nt protospacer, 1 = PAM-distal (5\' end of the spacer), increasing toward the PAM">',
  sprintf('##INFO=<ID=PS_REF,Number=1,Type=String,Description="Reference base on the PROTOSPACER strand (always %s here: the %s substrate). Equals REF on + sites, complements it on - sites">',
          spec$sub, spec$substrate),
  sprintf('##INFO=<ID=PS_ALT,Number=1,Type=String,Description="Edited base on the PROTOSPACER strand (always %s here). Equals ALT on + sites, complements it on - sites">',
          spec$edit),
  sprintf('##INFO=<ID=PS_STRAND,Number=1,Type=String,Description="Strand of the scored protospacer. On - the protospacer-strand %s>%s edit is written REF=%s ALT=%s, the same event seen from the reference strand">',
          spec$sub, spec$edit, COMPLEMENT[[spec$sub]], COMPLEMENT[[spec$edit]]),
  '##INFO=<ID=PS_SITE,Number=1,Type=String,Description="Scored protospacer interval, chrom:start-end, 1-based inclusive (spacer only, PAM excluded)">',
  '##INFO=<ID=PS_SOURCE,Number=1,Type=String,Description="Which finder located the site: seed (step 6, 6A_guidescan.py exact-seed anchor) or enumerate (step 7, GuideScan2, no seed requirement)">',
  '##INFO=<ID=PS_MM,Number=1,Type=Integer,Description="Spacer mismatches of the site alignment">',
  '##INFO=<ID=PS_RNAB,Number=1,Type=Integer,Description="RNA bulges in the site alignment">',
  '##INFO=<ID=PS_DNAB,Number=1,Type=Integer,Description="DNA bulges in the site alignment">',
  '##INFO=<ID=PS_CFD,Number=1,Type=Float,Description="Doench 2016 CFD score of the site alignment">',
  '##INFO=<ID=PS_ONTGT,Number=0,Type=Flag,Description="Site is the on-target (0 mismatch, 0 bulge) protospacer. 9_protospacer.R excludes these from its bar plots">',
  '##INFO=<ID=PS_PEAK,Number=1,Type=String,Description="Differential peak containing the site, chrom:start-end, 1-based inclusive">'
)

info_lines_possible <- function(spec) c(
  sprintf('##INFO=<ID=PS_OBS,Number=0,Type=Flag,Description="A %s>%s edit was actually observed at this available %s (this record also appears in the matching *_observed_%s_edits.vcf)">',
          spec$sub, spec$edit, spec$substrate, spec$tag),
  '##INFO=<ID=PS_AF,Number=1,Type=Float,Description="Mean FORMAT/AF across the 3 treatment replicates of the observed edit; absent unless PS_OBS">'
)

info_str <- function(.n, ...) {
  kv <- list(...)
  n <- .n
  if (n == 0) return(character(0))
  parts <- lapply(names(kv), function(k) {
    v <- rep_len(kv[[k]], n)
    if (is.logical(v)) ifelse(!is.na(v) & v, k, NA_character_)
    else ifelse(is.na(v), NA_character_, paste0(k, "=", v))
  })
  m <- matrix(unlist(parts), nrow = n)
  apply(m, 1, function(r) {
    r <- r[!is.na(r)]
    if (length(r) == 0) "." else paste(r, collapse = ";")
  })
}

fmt_num <- function(x, digits = 6) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "g", digits = digits))
}

region_label <- function(chrom, start0, end0) sprintf("%s:%d-%d", chrom, start0 + 1L, end0)

write_vcf <- function(path, meta, chrom_line, body) {
  con <- file(path, "wt")
  on.exit(close(con))
  writeLines(c(meta, chrom_line), con)
  if (length(body) > 0) writeLines(body, con)
}

summaries <- list()

for (cond in conditions) {
  cat(sprintf("\n=== %s ===\n", cond))
  tsv_path <- file.path(seed_dir, seed_files[[cond]])
  vcf_path <- file.path(vcf_dir, cond,
                        sprintf("%s_vs_control_relcontrol_filt_no-alt-chr.vcf", cond))
  samples <- tumor_samples[[cond]]
  stopifnot(file.exists(tsv_path), file.exists(vcf_path))

  # editor-independent: sites, spacer sequences and variant/site join, shared by both passes
  peaks <- read_sites_all(tsv_path, cond)
  sites <- peaks[peaks$guidescan_match, ]
  sites$on_target <- !is.na(sites$mismatches) & sites$mismatches == 0L &
    sites$rna_bulges == 0L & sites$dna_bulges == 0L
  sites$site_id <- seq_len(nrow(sites))
  if (sum(sites$on_target) > 1 ||
      (sum(sites$on_target) == 0 && !ALLOW_NO_ON_TARGET)) {
    stop(sprintf("%s: expected exactly 1 on-target site, found %d",
                 cond, sum(sites$on_target)))
  }
  cat(sprintf("  %d peaks, %d with a scored protospacer (%d seed, %d enumerate), %d unscored\n",
              nrow(peaks), nrow(sites), sum(sites$site_source == "seed"),
              sum(sites$site_source == "enumerate"), sum(!peaks$guidescan_match)))

  site_label <- region_label(sites$site_chrom, sites$site_start, sites$site_end)
  peak_label <- region_label(sites$peak_chrom, sites$peak_start, sites$peak_end)

  # site x position base matrix, bulges handled as in substrate_by_position
  spacer <- site_spacer_sequences(sites)
  chars <- matrix(NA_character_, nrow(sites), PROTOSPACER_LEN)
  for (i in seq_along(spacer)) {
    cs <- strsplit(spacer[i], "", fixed = TRUE)[[1]]
    k <- min(length(cs), PROTOSPACER_LEN)
    if (k > 0) chars[i, seq_len(k)] <- cs[seq_len(k)]
  }

  vr <- read_variants_raw(vcf_path, samples)
  variants <- vr$df

  sites_gr <- GRanges(seqnames = sites$site_chrom,
                      ranges = IRanges(start = sites$site_start + 1L, end = sites$site_end),
                      strand = sites$strand)
  variants_gr <- GRanges(seqnames = variants$chrom,
                         ranges = IRanges(start = variants$pos, width = 1L))
  # select = "first": one site per variant, so a variant in two overlapping
  # sites yields one record -- the same join, and the same count, as the plots.
  hit <- findOverlaps(variants_gr, sites_gr, select = "first")
  m <- which(!is.na(hit))
  si <- hit[m]                                  # row in `sites` for each match
  cat(sprintf("  %d/%d substitutions (%.1f%%) matched a scored protospacer\n",
              length(m), nrow(variants),
              if (nrow(variants) > 0) 100 * length(m) / nrow(variants) else 0))

  s_strand <- sites$strand[si]
  ps_pos <- ifelse(s_strand == "+",
                   variants$pos[m] - (sites$site_start[si] + 1L) + 1L,
                   sites$site_end[si] - variants$pos[m] + 1L)
  # NTS-oriented REF/ALT (complemented on '-' sites), not collapsed
  ps_ref <- ifelse(s_strand == "-", unname(COMPLEMENT[variants$ref[m]]), variants$ref[m])
  ps_alt <- ifelse(s_strand == "-", unname(COMPLEMENT[variants$alt[m]]), variants$alt[m])

  # per-editor passes
  for (ed in editors) {
    spec <- EDITOR_SPECS[[ed]]
    sig <- sprintf("%s>%s", spec$sub, spec$edit)
    cat(sprintf("\n  --- %s (protospacer-strand %s) ---\n", ed, sig))

    # OBSERVED: this editor's signature edits inside a scored protospacer
    keep <- which(ps_ref == spec$sub & ps_alt == spec$edit)
    outside <- sum(ps_pos[keep] < 1L | ps_pos[keep] > PROTOSPACER_LEN)
    if (outside > 0) {
      warning(sprintf("%s/%s: %d %s edit(s) outside protospacer positions 1-%d (bulged site longer than the spacer); they cannot appear in the possible file",
                      cond, ed, outside, sig, PROTOSPACER_LEN))
      keep <- keep[ps_pos[keep] >= 1L & ps_pos[keep] <= PROTOSPACER_LEN]
    }

    o <- list(
      site_id = sites$site_id[si][keep], ps_pos = ps_pos[keep],
      chrom = variants$chrom[m][keep], pos = variants$pos[m][keep],
      af = variants$af[m][keep], on_target = sites$on_target[si][keep]
    )
    ord <- order(match(o$chrom, fai$contig), o$pos)

    obs_info <- info_str(
      length(ord),
      PS_COND = cond, PS_EDITOR = ed, PS_POS = o$ps_pos[ord],
      PS_REF = spec$sub, PS_ALT = spec$edit,
      PS_STRAND = sites$strand[si][keep][ord],
      PS_SITE   = site_label[o$site_id[ord]],
      PS_SOURCE = sites$site_source[o$site_id[ord]],
      PS_MM     = sites$mismatches[o$site_id[ord]],
      PS_RNAB   = sites$rna_bulges[o$site_id[ord]],
      PS_DNAB   = sites$dna_bulges[o$site_id[ord]],
      PS_CFD    = fmt_num(sites$cfd[o$site_id[ord]]),
      PS_ONTGT  = o$on_target[ord],
      PS_PEAK   = peak_label[o$site_id[ord]]
    )

    vi <- m[keep][ord]     # rows of `variants` in output order
    obs_body <- paste(variants$chrom[vi], variants$pos[vi], variants$id[vi],
                      variants$ref[vi], variants$alt[vi], variants$qual[vi],
                      variants$filter[vi], obs_info, variants$fmt[vi],
                      vapply(variants$smp[vi], paste, character(1), collapse = "\t"),
                      sep = "\t")

    obs_path <- file.path(vcf_out_dir, sprintf("%s_observed_%s_edits.vcf", cond, spec$tag))
    write_vcf(obs_path,
              meta = c(vr$meta, info_lines_common(spec),
                       sprintf('##source_10_vcf="10_vcf.R observed %s protospacer edits, editor=%s, condition=%s"',
                               sig, ed, cond)),
              chrom_line = paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                                   "INFO", "FORMAT", samples), collapse = "\t"),
              body = obs_body)

    # POSSIBLE: every substrate base available to this editor
    idx <- which(chars == spec$sub, arr.ind = TRUE)   # NA compares NA, so ignored
    idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
    p_site <- idx[, 1]
    p_pos  <- idx[, 2]

    plus <- sites$strand[p_site] == "+"
    # protospacer position -> genomic coordinate and reference-strand allele
    p_gpos <- ifelse(plus, sites$site_start[p_site] + p_pos,
                     sites$site_end[p_site] - p_pos + 1L)
    p_ref <- ifelse(plus, spec$sub,  COMPLEMENT[[spec$sub]])
    p_alt <- ifelse(plus, spec$edit, COMPLEMENT[[spec$edit]])

    # flag the ones where an edit was actually seen, by (site, position)
    okey <- paste(o$site_id, o$ps_pos)
    pkey <- paste(sites$site_id[p_site], p_pos)
    oi <- match(pkey, okey)
    p_obs <- !is.na(oi)
    p_af  <- ifelse(p_obs, fmt_num(o$af[oi]), NA_character_)

    pord <- order(match(sites$site_chrom[p_site], fai$contig), p_gpos)

    pos_info <- info_str(
      length(pord),
      PS_COND = cond, PS_EDITOR = ed, PS_POS = p_pos[pord],
      PS_REF = spec$sub, PS_ALT = spec$edit,
      PS_STRAND = sites$strand[p_site][pord],
      PS_SITE   = site_label[p_site][pord],
      PS_SOURCE = sites$site_source[p_site][pord],
      PS_MM     = sites$mismatches[p_site][pord],
      PS_RNAB   = sites$rna_bulges[p_site][pord],
      PS_DNAB   = sites$dna_bulges[p_site][pord],
      PS_CFD    = fmt_num(sites$cfd[p_site][pord]),
      PS_ONTGT  = sites$on_target[p_site][pord],
      PS_PEAK   = peak_label[p_site][pord],
      PS_OBS    = p_obs[pord],
      PS_AF     = p_af[pord]
    )

    pos_body <- paste(sites$site_chrom[p_site][pord], p_gpos[pord], ".",
                      p_ref[pord], p_alt[pord], ".", ".", pos_info, sep = "\t")

    pos_path <- file.path(vcf_out_dir, sprintf("%s_possible_%s_edits.vcf", cond, spec$tag))
    write_vcf(pos_path,
              meta = c("##fileformat=VCFv4.2",
                       sprintf("##reference=%s", genome_fa),
                       sprintf('##contig=<ID=%s,length=%d>', fai$contig, fai$length),
                       info_lines_common(spec), info_lines_possible(spec),
                       sprintf('##source_10_vcf="10_vcf.R all available protospacer %ss (%s substrate), editor=%s, condition=%s"',
                               spec$substrate, ed, ed, cond)),
              chrom_line = paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),
                                 collapse = "\t"),
              body = pos_body)

    # CHECKS (fatal). 1: every observed edit is in the possible set (validates strand/position)
    stopifnot(all(okey %in% pkey))
    # 2: observed REF matches the reference base
    n_check <- min(length(vi), 2000L)
    if (n_check > 0) {
      chk <- sample(seq_along(vi), n_check)
      reg <- tempfile(fileext = ".txt")
      writeLines(sprintf("%s:%d-%d", variants$chrom[vi][chk],
                         variants$pos[vi][chk], variants$pos[vi][chk]), reg)
      raw <- system2(SAMTOOLS, c("faidx", shQuote(genome_fa), "-r", shQuote(reg)),
                     stdout = TRUE, stderr = FALSE)
      unlink(reg)
      got <- toupper(raw[!startsWith(raw, ">")])
      stopifnot(length(got) == n_check, all(got == variants$ref[vi][chk]))
    }
    # 3. Every possible-file position must lie inside its own site's interval.
    stopifnot(all(p_gpos > sites$site_start[p_site] & p_gpos <= sites$site_end[p_site]))

    n_obs_off <- sum(!o$on_target)
    n_pos_off <- sum(!sites$on_target[p_site])
    dup_pos <- sum(duplicated(paste(sites$site_chrom[p_site], p_gpos)))

    cat(sprintf("    observed %s : %d records (%d off-target, %d on-target)\n",
                sig, length(o$pos), n_obs_off, sum(o$on_target)))
    cat(sprintf("    possible %s : %d records (%d off-target, %d on-target; %d duplicate genomic position(s) from overlapping sites)\n",
                sig, length(p_gpos), n_pos_off, sum(sites$on_target[p_site]), dup_pos))
    cat(sprintf("    off-target rate = %.2f per 1,000 available %ss\n",
                1000 * n_obs_off / max(1L, n_pos_off), spec$substrate))

    summaries[[paste(cond, ed)]] <- data.frame(
      editor = ed, signature = sig, condition = cond, sites = nrow(sites),
      observed = length(o$pos), observed_offtarget = n_obs_off,
      possible = length(p_gpos), possible_offtarget = n_pos_off,
      per_1k_offtarget = 1000 * n_obs_off / max(1L, n_pos_off),
      stringsAsFactors = FALSE)
  }
}

summary_df <- do.call(rbind, summaries)
summary_df <- summary_df[order(match(summary_df$editor, editors),
                               match(summary_df$condition, conditions)), ]
write.csv(summary_df, file.path(vcf_out_dir, "vcf_record_counts.csv"), row.names = FALSE)

cat("\n=== record counts ===\n")
print(summary_df, row.names = FALSE)

# ABE vs CBE rate
if (length(editors) > 1 && all(c("ABE", "CBE") %in% editors)) {
  cat("\n=== signature vs background (off-target, per 1,000 available substrate bases) ===\n")
  for (cond in conditions) {
    a <- summary_df[summary_df$condition == cond & summary_df$editor == "ABE", ]
    b <- summary_df[summary_df$condition == cond & summary_df$editor == "CBE", ]
    if (nrow(a) == 1 && nrow(b) == 1) {
      cat(sprintf("  %-14s A>G %7.2f   C>T %6.2f   ratio %6.1fx\n",
                  cond, a$per_1k_offtarget, b$per_1k_offtarget,
                  a$per_1k_offtarget / max(b$per_1k_offtarget, .Machine$double.eps)))
    }
  }
}

cat(sprintf("\nWrote %s\n", vcf_out_dir))
