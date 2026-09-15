#!/usr/bin/env Rscript
# =============================================================================
# 9_protospacer.R
# Variant counts by position within the gRNA protospacer, colored by
# substitution type, for each condition's background-subtracted calls from
# 8_gatk_subtract.sh.
#
# Adapted from 260709_HEK293T_editors_revs/Scripts/8_strandedness.R for this
# project's design: ONE editor (ABE8e) across FOUR loci/sort conditions
# (HEK4, HEK4_noPalbo, HEK4_noSort, CD7) vs. a shared "control" (background,
# no locus-targeting sgRNA) baseline -- not that script's six-editor-vs-dCas9
# comparison. Same core logic otherwise:
#
# 6_guidescan_HEK4.sh / 6_guidescan_CD7.sh's ../results/09_guidescan/<condition>_
# binding_sites_seed.tsv already gives each peak's scored protospacer
# footprint (chrom/start/end) and strand; 8_gatk_subtract.sh's background-
# subtracted VCFs are NOT restricted to just those coordinates (unlike the
# source project's 7_gatk.sh, which calls with -L against the binding-site
# BED) -- they're still scoped to each condition's full differential-peaks
# region (see mutect2_call_variants.nf). So the "% matched a scored
# protospacer" logged below is expected to be well under 100%, and only
# matched variants are plotted. About half of matched sites are on the minus
# strand, where the spacer's 5' end is at the HIGHER genomic coordinate --
# naively using (variant_pos - site_start) would silently swap "near the 5'
# end" and "near the PAM" for half the sites. This script resolves that with
# a strand-correct 1-based position (1 = 5' end / PAM-distal, increasing
# toward the PAM) before counting.
#
# Substitution COLOR uses the same 6-class pyrimidine-collapsed scheme used by
# MutationalPatterns' plot_96_profile (not run here -- see 8_mutational_
# patterns.R's dependency notes; this script only needs GenomicRanges/
# ggplot2/dplyr/patchwork, already installed except patchwork). Indels
# (REF/ALT not both a single base) have no substitution type and are
# excluded, logged alongside the protospacer-match rate.
#
# That collapse merges each substitution with its reverse complement, so on its
# own it cannot say WHICH DNA STRAND was edited -- A>G and T>C land in one
# bucket. This script therefore runs in TWO MODES, selected by the STRANDED
# flag defined below (see that block for the full list of what each changes);
# both write the same set of filenames, aggregated mode into its own
# subdirectory. With STRANDED = TRUE the bar plots are MIRRORED about zero and
# the strand is carried by the bar's direction rather than by extra colors (12
# categorical fills would exceed what is distinguishable; direction is free and
# reads instantly):
#
#   UP    the protospacer-oriented substitution already has a pyrimidine REF,
#         i.e. it IS the collapsed label (e.g. T>C = adenine deaminated on the
#         TARGET strand, so the protospacer base was the paired T)
#   DOWN  its purine-REF complement, same color class (e.g. A>G = adenine
#         deaminated on the PROTOSPACER/displaced strand -- ABE8e's expected
#         substrate, since that strand is the ssDNA exposed by the R-loop)
#
# "Protospacer-oriented" means rewritten onto the strand whose sequence matches
# the gRNA spacer; VCF REF/ALT are always on the reference (+) strand, so a
# minus-strand site has both bases complemented (see ps_substitution()).
# peak_snp_summary.csv's type_1..type_20 store this strand-resolved form
# directly rather than the collapsed class.
#
# Reads {condition}_vs_control_relcontrol_filt_no-alt-chr.vcf (all
# substitution types, not signature-filtered) and does not recompute any
# calls -- only bins them by protospacer position and type.
#
# Each condition's guidescan TSV contains exactly one on-target site (the
# unique 0-mismatch/0-bulge hit -- see read_sites()); that site's own SNP(s)
# are excluded from all three bar plots below (they'd dwarf the off-target
# distribution the plots exist to show) and reported separately instead:
# on-target AF +/- SEM is written to on_target_af.csv and summarized as text
# in each panel's title.
# =============================================================================

suppressMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# =============================================================================
# PLOT_FONT -- Arial, by way of Helvetica.
#
# Arial itself cannot be used directly: it is not installed in becaskas.sif, and
# neither is any metric-compatible substitute (Liberation Sans, Arimo and Nimbus
# Sans are all absent), nor is there a host font directory to bind in -- the
# container ships only DejaVu Sans, Ubuntu and Source Code Pro. Checked with
# `fc-list` inside the container and on the login node.
#
# "Helvetica" is the standard way out: it is one of the 14 built-in PDF base
# fonts, so it needs no font file at all, and it is metrically identical to
# Arial -- which is exactly why PDF viewers and Illustrator substitute one for
# the other by default. The PDFs written below therefore carry Arial metrics and
# open as Arial/Helvetica in a typesetting tool. Same reasoning, and the same
# choice, as 13_venn_hek4_vs_nopalbo.R's own font comment.
#
# CAVEAT -- this applies to the PDFs ONLY. The PNGs are rasterised through
# ragg/cairo, which needs a real font FILE on disk and will silently fall back
# to the container's default sans (DejaVu Sans) for an unmatched family. So use
# the PDFs for anything that will be typeset or figure-assembled; treat the PNGs
# as previews whose letterforms are not Arial. To make the PNGs match too, an
# Arial-metric TTF (Arimo or Liberation Sans) has to be installed and registered
# with systemfonts -- both packages are present in the container, only the font
# file is missing.
#
# update_geom_defaults covers every geom_text/geom_label in this script at once;
# the theme_minimal() calls pass base_family explicitly, since base_family
# defaults back to "" on each call and would otherwise override theme_set().
PLOT_FONT <- "Helvetica"
theme_set(theme_minimal(base_family = PLOT_FONT))
update_geom_defaults("text", list(family = PLOT_FONT))

# Replaces Scripts/9_protospacer.R's hardcoded paths and per-DESIGN literals with
# pipeline-supplied config. Everything below this block is unchanged from that
# script apart from two path lookups (enumerate_files, peakbed_files), so the port
# stays diffable against the reference. Only its `condition` design is ported.
suppressPackageStartupMessages(library(optparse))

.opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character",
              help = "JSON: groups, tumor_samples, seed_files, enumerate_files, peakbed_files, palette, shared_groups, panel_ncol, panel_w, panel_h"),
  make_option("--stranded", type = "character", default = "TRUE",
              help = "TRUE = strand-resolved mirrored bars; FALSE = aggregated [default TRUE]"),
  make_option("--seed-dir", type = "character", default = "."),
  make_option("--vcf-dir",  type = "character", default = "."),
  make_option("--diff-dir", type = "character", default = "."),
  make_option("--out-dir",  type = "character", default = "."),
  make_option("--genome",   type = "character", help = "Indexed reference FASTA"),
  make_option("--samtools", type = "character", default = "samtools")
)))

.cfg <- jsonlite::fromJSON(.opt$config, simplifyVector = FALSE)

seed_dir  <- .opt$`seed-dir`
vcf_dir   <- .opt$`vcf-dir`
diff_dir  <- .opt$`diff-dir`
# Read straight from the reference rather than parsed out of any tool's
# match_sequence encoding -- see site_spacer_sequences() below.
genome_fa <- .opt$genome
SAMTOOLS  <- .opt$samtools

# TRUE = strand-resolved, bars mirrored about zero; FALSE = every substitution
# pooled with its reverse complement, single-sided. See the README.
STRANDED <- as.logical(.opt$stranded)
if (is.na(STRANDED)) stop("--stranded must be TRUE or FALSE (got: ", .opt$stranded, ")")
cat(sprintf("STRANDED = %s (%s)\n", STRANDED,
            if (STRANDED) "strand-resolved, mirrored bars" else "aggregated, non-stranded"))

# Aggregated mode writes to its own subdirectory so the two modes can coexist.
out_dir <- if (STRANDED) .opt$`out-dir` else file.path(.opt$`out-dir`, "aggregated")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Stack stage 2's sites alongside stage 1's; `site_source` marks each row's origin.
# The two are not homogeneous -- different mismatch ceilings, different ranking.
INCLUDE_ENUMERATE <- TRUE

SIGN_LEVELS <- if (STRANDED) c(-1L, 1L) else 1L
ZERO_LINE <- if (STRANDED) geom_hline(yintercept = 0, linewidth = 0.4, color = "black") else NULL
mode_ylim <- function(lim) if (STRANDED) c(-lim, lim) else c(0, lim)

# --- per-group config, from --config ---------------------------------------
conditions    <- unlist(.cfg$groups)
# treatment replicates only -- the merged VCF also carries control columns
tumor_samples <- lapply(.cfg$tumor_samples, unlist)
seed_files      <- unlist(.cfg$seed_files)        # stage 1 binding-site TSV per group
enumerate_files <- unlist(.cfg$enumerate_files)   # stage 2 TSV per group
peakbed_files   <- unlist(.cfg$peakbed_files)     # differential peaks BED per group
condition_palette <- unlist(.cfg$palette)
# groups sharing one y-axis; only same-gRNA groups are comparable
shared_groups <- unlist(.cfg$shared_groups)
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

# CD7 is a different locus from the three HEK4-locus conditions and has far
# fewer matched variants, so forcing it onto the same y-axis maximum as
# HEK4/HEK4_noPalbo/HEK4_noSort (directly comparable to each other -- same
# locus, same potency scale) flattens its panel. `totals_df` needs
# `condition`/`total` columns (position_totals / af_position_totals /
# altcount_position_totals below all have that shape).
#
# The bar plots are mirrored around zero (see strand_sign_of below), so this
# returns a SYMMETRIC half-range: the largest stack height in either
# direction, applied as ylim = c(-lim, +lim). Symmetric on purpose -- an
# asymmetric mirrored axis makes an up-bar and a down-bar of equal magnitude
# look different, which is the one comparison these plots exist to make.
panel_ylim <- function(totals_df, cond, shared_conditions = shared_groups) {
  pick <- if (cond %in% shared_conditions) {
    totals_df$condition %in% shared_conditions
  } else {
    totals_df$condition == cond
  }
  lim <- max(abs(totals_df$total[pick]), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) 1 else lim
}

# collapse a substitution and its complement onto the same pyrimidine-based
# class. This is the COLOR class only -- it deliberately merges each
# substitution with its reverse-complement (e.g. A>G with T>C), so on its own
# it cannot say which DNA strand was edited. strand_sign_of() below recovers
# that and encodes it as the bar's direction instead of a 7th-12th hue.
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

# The substitution rewritten in PROTOSPACER-STRAND orientation -- the strand
# whose sequence matches the gRNA spacer, i.e. the non-target strand that is
# displaced as ssDNA in the Cas9 R-loop. VCF REF/ALT are always given on the
# reference (+) strand, so a site whose protospacer lies on '-' has both bases
# complemented. NOT collapsed: keeping A>G distinct from T>C is the whole
# point.
#
# For ABE8e the two forms mean physically different things:
#   A>G  the deaminated adenine was on the protospacer (displaced) strand
#   T>C  the deaminated adenine was on the target strand (so the protospacer
#        base was the paired T, which becomes C)
ps_substitution <- function(ref, alt, site_strand) {
  neg <- !is.na(site_strand) & site_strand == "-"
  r <- ref
  a <- alt
  r[neg] <- unname(COMPLEMENT[ref[neg]])
  a[neg] <- unname(COMPLEMENT[alt[neg]])
  paste0(r, ">", a)
}

# Bar direction. All six collapsed classes carry a pyrimidine (C/T) REF, so a
# protospacer-oriented substitution that already has a C/T REF *is* its
# collapsed label and is drawn UP (+1); its purine-REF complement -- same
# color class, opposite strand -- is drawn DOWN (-1). Concretely for ABE8e:
# a T>C bar above the axis is a target-strand edit, and the A>G bar below it
# is the protospacer/displaced-strand edit that the collapse used to hide.
strand_sign_of <- function(ps_type) {
  ifelse(substr(ps_type, 1, 1) %in% c("C", "T"), 1L, -1L)
}

# per-sample FORMAT/AF across `samples`' own columns for one variant line,
# given that line's own FORMAT string (which field order can vary row-to-row
# -- e.g. phased vs unphased records drop PGT/PID/PS -- so AF's index is
# looked up per line rather than assumed fixed). NA for a replicate with no
# call ("./.") or a missing AF field. Returned as the raw per-replicate
# vector (not collapsed to a mean here) because the on-target report below
# needs a proper across-replicate SEM, not just a single averaged number.
sample_afs <- function(format_str, sample_fields) {
  af_idx <- match("AF", strsplit(format_str, ":", fixed = TRUE)[[1]])
  vapply(sample_fields, function(s) {
    # AF is Number=A (one value per ALT allele, comma-separated); this script
    # only ever uses af for atomized single-ALT SNP rows, but multi-allelic
    # indel rows (e.g. a spanning-deletion "*" placeholder ALT, AF="0.13,.")
    # still pass through here before being dropped by the is_snp filter, so
    # take the first allele's value to avoid an NA-coercion warning on ".".
    raw <- strsplit(s, ":", fixed = TRUE)[[1]][af_idx]
    if (is.na(raw) || raw == ".") return(NA_real_)
    as.numeric(strsplit(raw, ",", fixed = TRUE)[[1]][1])
  }, numeric(1))
}

# standard error of the mean, ignoring NAs (missing/uncalled replicates);
# NA if fewer than 2 non-missing values (SD, and therefore SEM, undefined).
sem <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / sqrt(length(x))
}

# total ALT-supporting reads across `samples`' own columns for one variant
# line -- a raw read count (FORMAT/AD's alt component, summed across the 3
# replicates), unlike mean_af() above which is Mutect2's per-sample
# statistical AF estimate. AD is Number=R (ref then alt, one value per
# allele); same first-allele-only handling as mean_af for multi-allelic rows.
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
  # per-replicate AF, one numeric vector per row -- `af` below is just its
  # mean (used everywhere except the on-target report, which needs the
  # per-replicate values for a proper SEM; see position_in_protospacer()).
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

  # Mutect2 emits the same SNV twice when one base falls in two different
  # assembled haplotype events, each carrying its own AD/AF estimate -- this
  # originates in the per-sample calls (results/07_variants/<cond>/<sample>/,
  # 13-19 duplicated chrom/pos/ref/alt per replicate here), NOT in
  # 8_gatk_subtract.sh: the duplicates are present before its `bcftools norm
  # --atomize` step and are unchanged by it. Left in place, one edit would be
  # counted twice by table() in the count plot and summed twice in the ALT-read
  # plot (the AF plot is unaffected, since it averages identical-position
  # values).
  #
  # Collapse to one record per (chrom, pos, ref, alt), keeping the HIGHEST-AF
  # one. bcftools' own `norm -d exact` would keep the FIRST instead, which
  # here is systematically the weaker call (the low-AF record shares its AD
  # with a neighbouring linked allele), so the choice is made explicitly here
  # rather than upstream. Original row order is preserved -- only losing rows
  # are dropped -- since downstream indexing into this frame is positional.
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

# -----------------------------------------------------------------------------
# Step 7's sites, normalized onto step 6's coordinate convention so the two can
# be stacked into one site table (see INCLUDE_ENUMERATE above).
#
# Two things have to be reconciled before they can share a table:
#   1. COORDINATES. 6A writes the SPACER only (always 20 bp here). The
#      enumerate driver writes spacer+PAM, so 21-25 bp once bulges shift it
#      (span = 20 - rna + dna + 3). Trimming the 3 nt PAM is strand-aware --
#      the PAM sits at the HIGH-coordinate end on '+' and the LOW end on '-'.
#      Without this the merged sites would contribute PAM positions and the
#      1..20 protospacer axis would be wrong for every enumerate row.
#   2. seed_hits. Absent from step 7's schema; set to 0, which is exactly what
#      it means -- these are by construction the peaks where Stage A found no
#      seed anchor.
#
# Returns NULL (with a warning) when the file is missing, so a condition whose
# step-7 job has not finished yet degrades to step-6-only rather than failing.
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

# Both site sources for one condition, stacked. `site_source` is carried
# through so any row in the outputs can be traced back to which tool found it.
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
      # step 7 only ever runs on peaks step 6 scored as "." , which the filter
      # above already dropped, so this cannot double-count a peak.
      df <- rbind(df[, intersect(names(df), names(enum))],
                   enum[, intersect(names(df), names(enum))])
    }
  }
  # The on-target site is the unique perfect (0 mismatch, 0 bulge) hit in
  # this condition's guidescan output -- 6_guidescan_HEK4.sh / 6_guidescan_CD7.sh
  # must be run with SUMMIT_WIDTH=0 (full-peak search) for this to actually
  # be found; a narrower summit window can anchor the seed but leave no room
  # for Stage B to score it, silently missing the true on-target site.
  df$on_target <- as.integer(df$mismatches) == 0 &
    as.integer(df$rna_bulges) == 0 & as.integer(df$dna_bulges) == 0
  df[, c("chrom", "start", "end", "strand", "on_target", "site_source")]
}

position_in_protospacer <- function(vcf_path, tsv_path, samples, cond = NULL) {
  variants <- read_variants(vcf_path, samples)
  sites <- read_sites(tsv_path, cond)

  n_on_target <- sum(sites$on_target)
  if (n_on_target != 1) {
    stop(sprintf("%s: expected exactly 1 on-target (0 mismatch, 0 bulge) site, found %d",
                  tsv_path, n_on_target))
  }

  # sites: 0-based half-open (BED-style, see 6A_guidescan.py) -> 1-based closed.
  # These intervals are the 20 bp SPACER only -- 6A_guidescan.py's start/end
  # span the protospacer without the PAM (the PAM is present only in that
  # TSV's match_sequence column, making it 23 nt against a 20 bp interval).
  # So `position` below runs 1..20 and PAM-position variants never match.
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
  # site_id = row index into this condition's own `sites` table -- a stable
  # per-peak identifier (one scored site = one differential peak) used below
  # to count DISTINCT edited peaks, not just raw variant counts. on_target
  # flags the one site that is the intended target rather than an off-target
  # candidate -- used below to keep it out of the off-target bar plots and
  # report it separately instead (see on_target_summary).
  # ps_type / strand_sign can only be resolved HERE, not in read_variants():
  # they need the matched site's strand, which is not known until the overlap
  # join above. snp_type stays the collapsed 6-class color; ps_type is the
  # strand-resolved substitution and strand_sign is its plotting direction.
  ps_type <- ps_substitution(variants$ref[matched], variants$alt[matched], site_strand)

  result <- data.frame(position = position, snp_type = variants$snp_type[matched],
                        ps_type = ps_type,
                        # aggregated mode pins the sign to +1, which turns every
                        # mirroring step downstream into a no-op
                        strand_sign = if (STRANDED) strand_sign_of(ps_type) else 1L,
                        site_strand = site_strand,
                        site_source = sites$site_source[hit[matched]],
                        site_id = hit[matched], af = variants$af[matched],
                        alt_count = variants$alt_count[matched],
                        chrom = variants$chrom[matched], genomic_pos = variants$pos[matched],
                        ref = variants$ref[matched], alt = variants$alt[matched],
                        on_target = sites$on_target[hit[matched]])
  result$af_list <- variants$af_list[matched]
  result
}

# hits_by_condition[[cond]] holds EVERY matched substitution, on-target
# included; all_counts/all_af/all_altcount below are built from the
# off-target subset only (on_target == FALSE) -- the on-target site is not
# graphed in the bar plots (it would dwarf and dilute the off-target
# distribution they're meant to show). It's reported separately instead:
# see on_target_summary / on_target_af.csv and on_target_label below.
hits_by_condition <- list()
all_counts <- do.call(rbind, lapply(conditions, function(cond) {
  cat(sprintf("%s:\n", cond))
  vcf_path <- file.path(vcf_dir, cond, sprintf("%s_vs_control_relcontrol_filt_no-alt-chr.vcf", cond))
  tsv_path <- file.path(seed_dir, seed_files[[cond]])
  hits <- position_in_protospacer(vcf_path, tsv_path, tumor_samples[[cond]], cond)
  hits_by_condition[[cond]] <<- hits
  offtarget <- hits[!hits$on_target, ]
  offtarget$snp_type <- factor(offtarget$snp_type, levels = snp_types)
  # Cross-tabulate by strand_sign as well, then re-apply the sign to the count
  # so pyrimidine-REF (protospacer-oriented) substitutions stack UP and their
  # purine-REF complements stack DOWN. table() zero-fills every
  # position x type x sign cell, so a class absent on one strand still holds
  # its slot in the stack.
  offtarget$strand_sign <- factor(offtarget$strand_sign, levels = SIGN_LEVELS)
  as.data.frame(table(position = offtarget$position, snp_type = offtarget$snp_type,
                       strand_sign = offtarget$strand_sign)) %>%
    mutate(condition = cond,
           position = as.integer(as.character(position)),
           strand_sign = as.integer(as.character(strand_sign)),
           count = Freq * strand_sign) %>%
    select(-Freq)
}))
all_counts$condition <- factor(all_counts$condition, levels = conditions)

write.csv(all_counts, file.path(out_dir, "variant_position_in_protospacer.csv"), row.names = FALSE)

# =============================================================================
# On-target report: mean AF +/- SEM across the condition's 3 treatment
# replicates for every SNP found inside the on-target protospacer window
# (the 20 bp spacer only -- 6A_guidescan.py writes spacer-only start/end
# coordinates; the PAM appears in that TSV's match_sequence column but is
# NOT part of the interval, so a SNP in the PAM is outside this window and
# is not reported here).
# There can be more than one SNP in the window -- e.g. CD7 has a second,
# lower-AF edited adenine plus two near-zero-AF single-replicate calls that
# read as background noise. Written to its own CSV, kept out of the bar plots
# above, and summarized as text in each panel's title below.
# =============================================================================
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
    # Per-replicate AFs spelled out as "sample=value;sample=value;...". mean_af
    # and sem_af alone cannot say WHICH replicate carried a call, which matters
    # for the single-replicate rows (n_replicates = 1) that read as background
    # noise rather than real edits. `vals` is in tumor_samples[[cond]] order --
    # read_variants() builds it by indexing the VCF's sample columns in exactly
    # that order. This column is what 10_ontarget_af.py used to provide before
    # it was retired as redundant with this block.
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

# One label per condition for the plot titles below -- the on-target window
# can contain more than one SNP (see comment above), so the title uses the
# single highest-AF one (the primary intended edit); the full set, including
# any lower-confidence extras, is in on_target_af.csv.
on_target_label <- setNames(vapply(conditions, function(cond) {
  ot <- on_target_summary[on_target_summary$condition == cond, ]
  ot <- ot[!is.na(ot$mean_af), ]
  if (nrow(ot) == 0) return(sprintf("%s (on-target: not found)", cond))
  best <- ot[which.max(ot$mean_af), ]
  sprintf("%s (on-target AF %.2f\u00b1%.2f)", cond, best$mean_af, best$sem_af)
}, character(1)), conditions)

# Shared x range across all 4 panels. y range is shared among the three
# HEK4-locus conditions (directly comparable, same editor/locus) but CD7 --
# a different locus with far fewer matched variants -- gets its own, via
# panel_ylim() above (unlike the source script's 3 axis-groups for 6 editors
# of very different potency, this only splits out the one different locus).
# Totals are taken PER DIRECTION (up-stack vs down-stack) and not as a single
# signed sum -- summing signed counts would let an up-bar and a down-bar cancel
# and report a stack height of ~0 for a position that actually has tall bars in
# both directions. panel_ylim() then takes the largest absolute height.
signed_totals <- function(df, value_col) {
  df %>%
    mutate(dir = ifelse(.data[[value_col]] < 0, "down", "up")) %>%
    group_by(condition, position, dir) %>%
    summarise(total = sum(.data[[value_col]]), .groups = "drop")
}

# Per-class legend label. When STRANDED, show the up/down split so the strand
# asymmetry is legible as numbers and not only as bar geometry; when
# aggregated, one pooled total per class (a "/ dn 0" suffix would be noise).
signed_type_labels <- function(df, value_col, fmt = "%s") {
  v <- df[[value_col]]
  agg <- df %>%
    mutate(.v = v) %>%
    group_by(snp_type) %>%
    summarise(up = sum(.v[.v > 0]), down = -sum(.v[.v < 0]), .groups = "drop")
  if (STRANDED) {
    setNames(sprintf(paste0("%s (up ", fmt, " / dn ", fmt, ")"),
                      agg$snp_type, agg$up, agg$down), agg$snp_type)
  } else {
    setNames(sprintf(paste0("%s (", fmt, ")"), agg$snp_type, agg$up), agg$snp_type)
  }
}

position_totals <- signed_totals(all_counts, "count")
x_range <- range(all_counts$position)

# Kept to 3 short lines: a taller subtitle overruns the top row's panels.
MIRROR_SUBTITLE <- if (STRANDED) paste(
  "Off-target sites only; mirrored by edited strand, oriented to the protospacer (gRNA-matching) strand.",
  "UP = pyrimidine-REF form (T>C = adenine edited on the target strand).  DOWN = purine-REF complement",
  "(A>G = adenine edited on the protospacer/displaced strand, ABE8e's expected substrate).  Colors unchanged.",
  sep = "\n") else paste(
  "Off-target sites only. AGGREGATED (non-stranded) view: each substitution is pooled with its reverse",
  "complement into the 6 pyrimidine-collapsed classes, so A>G and T>C are counted together in T>C and the",
  "edited strand is NOT distinguished. See the strand-resolved version one directory up.",
  sep = "\n")

# Suffix appended to each plot's main title so the two modes are never confused
# when the PDFs are viewed side by side.
MODE_TAG <- if (STRANDED) ", mirrored by edited strand" else " (aggregated, non-stranded)"

# Axis titles are repeated in all 4 panels, so they must stay short or they
# collide with the neighbouring panel (and with the legend) -- the strand
# meaning lives in the subtitle above, not in the axis title.
X_LAB <- "Position in protospacer (1 = PAM-distal)"

# One ggplot per condition (not facet_wrap) so each panel's legend shows that
# condition's OWN per-type totals, not counts pooled across all 4 conditions.
plots <- lapply(conditions, function(cond) {
  df <- all_counts[all_counts$condition == cond, ]
  type_labels <- signed_type_labels(df, "count", "%d")
  y_lim <- panel_ylim(position_totals, cond)

  ggplot(df, aes(x = position, y = count, fill = snp_type)) +
    geom_col(color = "black", linewidth = 0.2, width = 0.8, alpha = 0.9) +
    ZERO_LINE +
    scale_fill_manual(values = snp_type_palette, labels = type_labels, name = "SNP type", drop = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = function(v) format(abs(v), big.mark = ",")) +
    coord_cartesian(xlim = x_range, ylim = mode_ylim(y_lim)) +
    labs(x = X_LAB, y = "Number of variants", title = on_target_label[[cond]]) +
    theme_minimal(base_size = 11, base_family = PLOT_FONT) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", color = condition_palette[[cond]]))
})

p <- wrap_plots(plots, ncol = panel_ncol) +
  plot_annotation(title = paste0("Variant counts by position in the gRNA protospacer, colored by SNP type", MODE_TAG),
                   subtitle = MIRROR_SUBTITLE,
                   theme = theme(text = element_text(family = PLOT_FONT), plot.title = element_text(face = "bold")))

ggsave(file.path(out_dir, "variant_position_in_protospacer.pdf"), p, width = panel_w, height = panel_h)
ggsave(file.path(out_dir, "variant_position_in_protospacer.png"), p, width = panel_w, height = panel_h, dpi = 200)

# =============================================================================
# Second version of the plot above: same x-axis (protospacer position) and
# same color-by-SNP-type stacking, but bar height is each type's AVERAGE
# allele fraction (FORMAT/AF, mean across that condition's own 3 treatment
# replicates -- see `tumor_samples` above; the control columns are excluded)
# instead of a raw variant count. Position/type combinations with no matched
# variant get an explicit AF of 0 so every panel still stacks across all 6
# types (mirroring how table() zero-fills the count version) -- but the
# per-type averages shown in each legend are computed only from actually-
# matched variants, not diluted by those structural zeros.
# =============================================================================
# Grid is now position x type x strand_sign, so a class present on only one
# strand still holds its slot in both stacks. AF is averaged within each
# (position, type, sign) cell and then signed, so the down-stack carries the
# same magnitudes it would have had before the split.
all_af <- do.call(rbind, lapply(conditions, function(cond) {
  hits <- hits_by_condition[[cond]]
  hits <- hits[!hits$on_target, ]
  hits$snp_type <- factor(hits$snp_type, levels = snp_types)
  if (nrow(hits) == 0) return(NULL)
  positions <- sort(unique(hits$position))
  grid <- expand.grid(position = positions,
                       snp_type = factor(snp_types, levels = snp_types),
                       strand_sign = SIGN_LEVELS)
  af_by_group <- hits %>% group_by(position, snp_type, strand_sign) %>%
    summarise(af = mean(af), .groups = "drop")
  grid %>%
    left_join(af_by_group, by = c("position", "snp_type", "strand_sign")) %>%
    mutate(af = ifelse(is.na(af), 0, af) * strand_sign, condition = cond)
}))
all_af$condition <- factor(all_af$condition, levels = conditions)

write.csv(all_af, file.path(out_dir, "variant_af_in_protospacer.csv"), row.names = FALSE)

af_position_totals <- signed_totals(all_af, "af")

af_plots <- lapply(conditions, function(cond) {
  df <- all_af[all_af$condition == cond, ]
  hits <- hits_by_condition[[cond]]
  hits <- hits[!hits$on_target, ]
  # Legend shows each class's mean AF per strand direction, from matched
  # variants only (not diluted by the structural zeros in `grid`).
  type_means <- hits %>% group_by(snp_type, strand_sign) %>%
    summarise(m = mean(af), .groups = "drop")
  lab <- vapply(snp_types, function(ty) {
    up <- type_means$m[type_means$snp_type == ty & type_means$strand_sign == 1L]
    dn <- type_means$m[type_means$snp_type == ty & type_means$strand_sign == -1L]
    f <- function(x) if (length(x) == 0 || is.na(x)) "-" else sprintf("%.2f", x)
    sprintf("%s (up %s / dn %s)", ty, f(up), f(dn))
  }, character(1))
  type_label_text <- setNames(lab, snp_types)
  af_y_lim <- panel_ylim(af_position_totals, cond)

  ggplot(df, aes(x = position, y = af, fill = snp_type)) +
    geom_col(color = "black", linewidth = 0.2, width = 0.8, alpha = 0.9) +
    ZERO_LINE +
    scale_fill_manual(values = snp_type_palette, labels = type_label_text, name = "SNP type", drop = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = function(v) sprintf("%.2f", abs(v))) +
    coord_cartesian(xlim = x_range, ylim = mode_ylim(af_y_lim)) +
    labs(x = X_LAB, y = "Average allele fraction", title = on_target_label[[cond]]) +
    theme_minimal(base_size = 11, base_family = PLOT_FONT) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", color = condition_palette[[cond]]))
})

p_af <- wrap_plots(af_plots, ncol = panel_ncol) +
  plot_annotation(title = paste0("Average allele fraction by position in the gRNA protospacer, colored by SNP type", MODE_TAG),
                   subtitle = MIRROR_SUBTITLE,
                   theme = theme(text = element_text(family = PLOT_FONT), plot.title = element_text(face = "bold")))

ggsave(file.path(out_dir, "variant_af_in_protospacer.pdf"), p_af, width = panel_w, height = panel_h)
ggsave(file.path(out_dir, "variant_af_in_protospacer.png"), p_af, width = panel_w, height = panel_h, dpi = 200)

# =============================================================================
# Third version of the same plot: bar height is now the TOTAL number of
# mutated alleles observed -- FORMAT/AD's alt component (raw supporting read
# count), summed across the condition's 3 treatment replicates per variant,
# then summed again across every distinct off-target site sharing a given
# position+type (unlike the AF plot, which averages at both steps).
# =============================================================================
all_altcount <- do.call(rbind, lapply(conditions, function(cond) {
  hits <- hits_by_condition[[cond]]
  hits <- hits[!hits$on_target, ]
  hits$snp_type <- factor(hits$snp_type, levels = snp_types)
  if (nrow(hits) == 0) return(NULL)
  positions <- sort(unique(hits$position))
  grid <- expand.grid(position = positions,
                       snp_type = factor(snp_types, levels = snp_types),
                       strand_sign = SIGN_LEVELS)
  altcount_by_group <- hits %>% group_by(position, snp_type, strand_sign) %>%
    summarise(alt_count = sum(alt_count), .groups = "drop")
  grid %>%
    left_join(altcount_by_group, by = c("position", "snp_type", "strand_sign")) %>%
    mutate(alt_count = ifelse(is.na(alt_count), 0, alt_count) * strand_sign, condition = cond)
}))
all_altcount$condition <- factor(all_altcount$condition, levels = conditions)

write.csv(all_altcount, file.path(out_dir, "variant_altcount_in_protospacer.csv"), row.names = FALSE)

altcount_position_totals <- signed_totals(all_altcount, "alt_count")

altcount_plots <- lapply(conditions, function(cond) {
  df <- all_altcount[all_altcount$condition == cond, ]
  type_labels <- signed_type_labels(df, "alt_count", "%d")
  altcount_y_lim <- panel_ylim(altcount_position_totals, cond)

  ggplot(df, aes(x = position, y = alt_count, fill = snp_type)) +
    geom_col(color = "black", linewidth = 0.2, width = 0.8, alpha = 0.9) +
    ZERO_LINE +
    scale_fill_manual(values = snp_type_palette, labels = type_labels, name = "SNP type", drop = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = function(v) format(abs(v), big.mark = ",")) +
    coord_cartesian(xlim = x_range, ylim = mode_ylim(altcount_y_lim)) +
    labs(x = X_LAB, y = "Total ALT reads (3 reps)", title = on_target_label[[cond]]) +
    theme_minimal(base_size = 11, base_family = PLOT_FONT) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", color = condition_palette[[cond]]))
})

p_altcount <- wrap_plots(altcount_plots, ncol = panel_ncol) +
  plot_annotation(title = paste0("Total mutated alleles (ALT reads summed across 3 replicates) by position in the gRNA protospacer", MODE_TAG),
                   subtitle = MIRROR_SUBTITLE,
                   theme = theme(text = element_text(family = PLOT_FONT), plot.title = element_text(face = "bold")))

ggsave(file.path(out_dir, "variant_altcount_in_protospacer.pdf"), p_altcount, width = panel_w, height = panel_h)
ggsave(file.path(out_dir, "variant_altcount_in_protospacer.png"), p_altcount, width = panel_w, height = panel_h, dpi = 200)

# =============================================================================
# Separate graphs: (1) raw number of differential peaks per condition, and
# (2) differential peaks per signature-matched (T>C, ABE8e's expected class)
# edit found -- lower means editing turns up more often per candidate site.
# Peak counts come straight from each condition's differential_peaks.bed
# (wc -l), not a separately-tracked counts file (this project has no
# equivalent of the source project's differential_peaks_counts.csv).
# =============================================================================
# diff_dir comes from --diff-dir (see the config block at the top)
peak_counts <- data.frame(
  condition = conditions,
  n_peaks = vapply(conditions, function(cond) {
    length(readLines(file.path(diff_dir, peakbed_files[[cond]])))
  }, integer(1))
)
peak_counts$condition <- factor(peak_counts$condition, levels = conditions)

write.csv(peak_counts, file.path(out_dir, "differential_peaks_per_condition.csv"), row.names = FALSE)

p_peak_counts <- ggplot(peak_counts, aes(x = condition, y = n_peaks, fill = condition)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7, alpha = 0.9) +
  geom_text(aes(label = format(n_peaks, big.mark = ",")), vjust = -0.3, size = 3) +
  scale_fill_manual(values = condition_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Differential peaks", title = "Differential peaks called, by condition") +
  theme_minimal(base_size = 12, base_family = PLOT_FONT) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "differential_peaks_per_condition.pdf"), p_peak_counts, width = 7, height = 5.5)
ggsave(file.path(out_dir, "differential_peaks_per_condition.png"), p_peak_counts, width = 7, height = 5.5, dpi = 200)

# The signature is now the STRAND-RESOLVED A>G -- an adenine deaminated on the
# protospacer (gRNA-matching, displaced) strand, which is the strand exposed as
# ssDNA in the R-loop and so ABE8e's expected substrate. The old collapsed
# "T>C" class pooled that with target-strand edits; n_collapsed_TC keeps the
# pooled figure alongside for continuity with earlier versions of this plot.
signature_edit_counts <- do.call(rbind, lapply(conditions, function(cond) {
  hits <- hits_by_condition[[cond]]
  hits <- hits[!hits$on_target, ]
  # Both figures are always reported; only which one drives the ratio changes.
  data.frame(condition = cond,
             signature = if (STRANDED) "A>G (protospacer strand)" else "T>C (aggregated, both strands)",
             n_signature_edits = if (STRANDED) sum(hits$ps_type == "A>G") else sum(hits$snp_type == "T>C"),
             n_strand_resolved_AG = sum(hits$ps_type == "A>G"),
             n_collapsed_TC = sum(hits$snp_type == "T>C"))
}))

peaks_per_edit <- merge(peak_counts, signature_edit_counts, by = "condition")
peaks_per_edit$condition <- factor(peaks_per_edit$condition, levels = conditions)
peaks_per_edit <- peaks_per_edit[order(peaks_per_edit$condition), ]
peaks_per_edit$peaks_per_edit <- ifelse(peaks_per_edit$n_signature_edits > 0,
                                         peaks_per_edit$n_peaks / peaks_per_edit$n_signature_edits, NA)

write.csv(peaks_per_edit, file.path(out_dir, "peaks_per_edit.csv"), row.names = FALSE)

p_ratio <- ggplot(peaks_per_edit, aes(x = condition, y = peaks_per_edit, fill = condition)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%s / %s = %s", format(n_peaks, big.mark = ","), n_signature_edits,
                                  ifelse(is.na(peaks_per_edit), "NA", sprintf("%.1f", peaks_per_edit)))),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = condition_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL,
       y = if (STRANDED) "Differential peaks per A>G edit (protospacer strand)"
            else "Differential peaks per T>C signature-matched edit",
       title = "Differential peaks per edit, by condition",
       subtitle = if (STRANDED) "Signature is the strand-resolved A>G (adenine deaminated on the displaced strand)"
                   else "Aggregated (non-stranded): T>C pools A>G and T>C, so the edited strand is not distinguished") +
  theme_minimal(base_size = 12, base_family = PLOT_FONT) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "peaks_per_edit.pdf"), p_ratio, width = 7, height = 5.5)
ggsave(file.path(out_dir, "peaks_per_edit.png"), p_ratio, width = 7, height = 5.5, dpi = 200)

# =============================================================================
# Per-peak wide summary table: ONE ROW PER DIFFERENTIAL PEAK (every peak, not
# just the ones with a hit), recording whether guidescan found a binding site
# there and, if so, the allele fraction and substitution type of any SNP at
# each protospacer position 1-20.
#
# This deliberately does NOT reuse hits_by_condition: read_sites() drops the
# unscored ('.') rows, so that path cannot represent a peak with no guidescan
# match -- which is half of what this table is for. read_sites_all() below
# keeps every row of the guidescan TSV instead. Each TSV row already IS one
# differential peak (6A_guidescan.py emits one row per input BED interval, in
# BED order, NA-filled when nothing scored), so peak_chrom/peak_start/
# peak_end is present even on unmatched rows and needs no join back to the
# original .bed.
#
# Columns:
#   condition, peak_chrom/peak_start/peak_end   the differential peak
#   guidescan_match                             TRUE if a site was scored here
#   seed_hits                                   Stage A anchors found (0 = no
#                                                 seed match anywhere in peak)
#   site_chrom/site_start/site_end/strand       the scored site (NA if none)
#   mismatches/rna_bulges/dna_bulges/cfd        its alignment stats (NA if none)
#   on_target                                   TRUE for the one perfect hit
#   n_snps                                      matched SNPs in this peak
#   af_1 .. af_20                               allele fraction at each
#                                                 strand-correct protospacer
#                                                 position (NA = no SNP there)
#   type_1 .. type_20                           that SNP's STRAND-RESOLVED
#                                                 substitution, oriented to the
#                                                 protospacer strand (A>G and
#                                                 T>C stay distinct -- NOT the
#                                                 6-class collapse used for
#                                                 plot fill color)
#
# Positions are the same strand-correct 1-based coordinate the plots use
# (1 = 5'/PAM-distal, increasing toward the PAM), so af_1..af_20 line up
# across plus- and minus-strand sites. AF is the mean across that condition's
# own 3 treatment replicates (same `af` the AF plot uses).
#
# On-target rows are INCLUDED here (unlike the bar plots, which exclude them)
# since this is a per-peak inventory rather than an off-target distribution --
# filter on on_target == FALSE to match the plots.
#
# Tie rule (now a fallback, not the primary defence): two SNPs could share one
# protospacer position only if a single genomic base carried two records.
# read_variants() already collapses those to the highest-AF record, so in
# practice no tie reaches here -- every site is exactly 20 bp (6A_guidescan.py's
# bulges always occur in equal rna/dna pairs, so 20 + dna - rna cancels), which
# means positions are always 1..20 and one base maps to one column. The rule is
# kept as a guard: if a tie ever did occur the higher-AF allele wins, and
# n_snps counting every matched SNP means n_snps exceeding the number of
# non-NA af_* columns flags any peak where that happened.
# =============================================================================
PROTOSPACER_LEN <- 20L

# Protospacer-strand sequence of each site, read from the REFERENCE over that
# site's own interval -- deliberately not parsed out of match_sequence.
#
# Why: the substrate denominator has to use the same notion of "position p" the
# plots use, and the plots are bulge-naive (position = variant_pos - site_start
# + 1, straight arithmetic on the interval). Reading the genome over the same
# interval reproduces that exactly. Parsing match_sequence would not: step 6
# and step 7 encode bulges differently, and for a DNA-bulge row guidescan.py
# itself notes the inserted base cannot be identified ("we cannot tell which
# lowercase base is the insertion"). Going to the reference sidesteps both
# problems and makes seed and enumerate sites strictly comparable.
#
# One batched `samtools faidx -r` per condition. Minus-strand sites are
# reverse-complemented so index 1 is always the PAM-distal end, matching the
# plots' position convention.
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

# Every row of the guidescan TSV, unscored ones included -- see section header.
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
    # 6A_guidescan.py writes match_sequence in PROTOSPACER orientation (it
    # slices the same strand-specific `target` string the alignment used), so
    # an 'A' here is a protospacer-strand adenine -- i.e. actual ABE8e
    # substrate -- with no further complementing needed. Verified against the
    # on-target rows, whose match_sequence reproduces the input gRNA exactly.
    # 23 nt: 20 nt spacer + 3 nt PAM.
    match_sequence  = ifelse(df$start == ".", NA_character_, df$match_sequence),
    site_source     = ifelse(df$start == ".", NA_character_, "seed"),
    stringsAsFactors = FALSE
  )

  # Fill step 6's unscored ('.') rows in from step 7 where it found a site.
  # FILLED, not appended -- every no-seed peak already has a row here, so
  # appending would duplicate the peak.
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
        # tie rule: keep the higher-AF allele (see section header)
        if (is.na(af_mat[row, pos]) ||
            (!is.na(variants$af[v]) && variants$af[v] > af_mat[row, pos])) {
          af_mat[row, pos] <- variants$af[v]
          # STRANDED: the protospacer-oriented substitution, keeping A>G
          # (adenine edited on the displaced strand) distinct from T>C (adenine
          # edited on the target strand). Aggregated: the 6-class collapse,
          # which merges that pair and so does not identify the strand.
          type_mat[row, pos] <- if (STRANDED) {
            ps_substitution(variants$ref[v], variants$alt[v], peaks$strand[row])
          } else {
            variants$snp_type[v]
          }
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

# =============================================================================
# Summary bar graph: differential peaks -> GuideScan matches -> edited R-loops,
# per condition, with the last two bars STACKED by which finder located the
# site. Derived entirely from peak_snp_summary above (one row per peak), so the
# bars cannot disagree with that table.
#
# Four bar positions per condition, left to right:
#   1. Differential peaks   every peak. NOT stacked -- an unmatched peak has no
#                           site and therefore no source, so the split is
#                           undefined here rather than zero.
#   2. GuideScan matches    stacked: grep (step 6, exact-seed-anchored) +
#                           enumerate (step 7, no seed requirement).
#   3. R-loops with edits   the same split, restricted to peaks carrying >=1 SNP
#                           inside the scored 20 bp protospacer.
#   4. Edited positions     the same split again, but counting DISTINCT EDITED
#                           GENOMIC POSITIONS rather than peaks.
#
# Bars 1-3 are nested peak counts, each a subset of the one before it. BAR 4 IS
# NOT: it changes the unit of the count from peaks to edited bases, so it
# necessarily EXCEEDS bar 3 by however many peaks carry more than one edit
# (795 vs 476 in HEK4, 1088 vs 712 in HEK4_noSort). The funnel therefore stops
# narrowing at bar 3 -- that widening is the point of the bar, not an error, and
# the y-axis title and subtitle both name the unit change so the reader is not
# invited to read bar 4 as a further subset.
#
# Bar 4 is counted as the number of non-NA af_1..af_20 cells per peak, i.e. one
# per protospacer position carrying an edit. That is a distinct GENOMIC position
# by construction (position p of a given site is one genomic base), and it comes
# from the same peak_snp_summary table as bars 1-3, so bar 4 >= bar 3 is
# guaranteed rather than hoped for -- the stopifnot below asserts it. It also
# reconciles exactly with the "matched a scored protospacer" counts logged by
# position_in_protospacer() further up (795 / 596 / 1088 / 63).
#
# The two segments are DISJOINT and sum to their bar's total: step 7 only ever
# ran on the peaks step 6 emitted as unscored ('.'), so no peak can be counted
# in both. `site_source` (set by read_sites_all) carries the split; the
# stopifnot below enforces the sum so a future change to either reader cannot
# silently desync the stack from the total.
#
# Why draw the split at all: the enumerate share is not a rounding detail (~3%
# of matches in HEK4, ~40% in HEK4_noSort), and the two populations do not
# behave alike -- enumerate sites are markedly ENRICHED for edits where the
# no-seed set is small (HEK4 51% vs 15% of grep sites; CD7 43% vs 17%) and
# markedly DEPLETED where it is large (HEK4_noSort 4.4% vs 10.7%). Pooling them
# into one bar, as this plot did previously, hides both facts -- and the second
# one matters, because a condition whose matches are mostly enumerate sites has
# its substrate denominator inflated by sites that carry no edits.
#
# With INCLUDE_ENUMERATE = FALSE the enumerate segments are 0 and this degrades
# to the previous three-bar figure.
#
# ENCODING. Metric = hue, in the reference categorical order (blue, orange,
# aqua, yellow = slots 1-4). Source = a two-step ordinal ramp WITHIN that
# metric's own hue, so the two segments read as one family rather than as two
# unrelated series; in every stacked bar the DARKER step is grep and the
# LIGHTER step is enumerate. Each pair was checked with the palette validator's
# --ordinal gate (monotone lightness, adjacent delta-L >= 0.06, light end >=
# 2:1 against the light surface): orange #eb6834/#f29a72 and aqua
# #1baf7a/#4dc294 pass, and yellow uses slot 4's own documented light and
# dark-mode steps, #eda100 over #c98500 -- a lightened yellow cannot be used
# here, since #eda100 is already only 2.11:1 against white and anything lighter
# fails the 2:1 light-end floor. Segments are separated by a surface-colored gap
# rather than the black hairline used elsewhere in this script -- with two steps
# of ONE hue touching, a dark outline reads as a third boundary rather than as a
# break between fills.
#
# Palette caveat: slots 1-4 pass every gate on the ADJACENT pairlist that bar
# and stack charts use (worst adjacent CVD delta-E 9.1, normal-vision 22.9), and
# the bar order puts yellow next to aqua rather than next to orange. Under the
# all-pairs pairlist the orange<->yellow pair would fail the normal-vision floor
# (13.7 < 15), so if this ever becomes a scatter/dot form rather than bars, the
# 4th slot has to be reconsidered. The mandatory direct labels below are also
# the documented relief for aqua (2.74:1) and yellow (2.11:1) sitting under 3:1.
#
# Linear y on purpose despite the ~51x spread (CD7 173 vs HEK4_noSort 8807):
# that spread is a real result, and every bar carries a direct total label so
# the short CD7 bars stay quantitative. Those labels also discharge the
# contrast relief rule for the aqua fills, which sit below 3:1 on white.
# =============================================================================
funnel_levels <- c("Differential peaks", "GuideScan matches", "R-loops with edits",
                    "Edited genomic positions")

# fill key = metric x source, listed in stack order (grep before enumerate, so
# grep ends up at the BOTTOM of each stack -- see position_stack(reverse=TRUE))
funnel_fills <- c("Differential peaks",
                  "GuideScan match: grep (step 6)",
                  "GuideScan match: enumerate (step 7)",
                  "Edited R-loop: grep (step 6)",
                  "Edited R-loop: enumerate (step 7)",
                  "Edited position: grep (step 6)",
                  "Edited position: enumerate (step 7)")
funnel_palette <- setNames(
  c("#2a78d6",               # blue        -- differential peaks (unstacked)
    "#eb6834", "#f29a72",    # orange ramp -- guidescan matches
    "#1baf7a", "#4dc294",    # aqua ramp   -- edited R-loops
    "#c98500", "#eda100"),   # yellow ramp -- edited genomic positions
  funnel_fills)

af_cols <- sprintf("af_%d", seq_len(PROTOSPACER_LEN))

funnel_counts <- do.call(rbind, lapply(conditions, function(cond) {
  d <- peak_snp_summary[peak_snp_summary$condition == cond, ]
  # %in% rather than == so the NA site_source on unmatched rows yields FALSE
  # instead of NA (which would poison the sums below).
  from_grep <- d$site_source %in% "seed"
  from_enum <- d$site_source %in% "enumerate"
  edited <- d$n_snps > 0
  # one per protospacer position carrying an edit = one distinct genomic base
  n_pos <- rowSums(!is.na(d[, af_cols]))
  stopifnot(sum(from_grep | from_enum) == sum(d$guidescan_match),
            # bar 4 counts bases, bar 3 counts peaks, so it can only be larger
            sum(n_pos) >= sum(edited))
  data.frame(
    condition = cond,
    metric = factor(funnel_levels[c(1, 2, 2, 3, 3, 4, 4)], levels = funnel_levels),
    site_source = c(NA, rep(c("grep", "enumerate"), 3)),
    fill = factor(funnel_fills, levels = funnel_fills),
    n = c(nrow(d),
          sum(from_grep), sum(from_enum),
          sum(edited & from_grep), sum(edited & from_enum),
          sum(n_pos[from_grep]), sum(n_pos[from_enum])),
    stringsAsFactors = FALSE)
}))
funnel_counts$condition <- factor(funnel_counts$condition, levels = conditions)
# each row as a % of that condition's own peak total (its first row)
funnel_counts$pct_of_peaks <- round(
  100 * funnel_counts$n / ave(funnel_counts$n, funnel_counts$condition,
                               FUN = function(v) v[1]), 1)

write.csv(funnel_counts, file.path(out_dir, "peak_guidescan_edit_counts.csv"), row.names = FALSE)

# ggplot cannot dodge and stack in a single geom, so the dodge is done by hand:
# each metric sits at a fixed offset from its condition's integer x position,
# and the segments that share an x stack. reverse = TRUE puts the first fill
# level (grep) at the bottom of the stack rather than the top.
FUNNEL_DODGE <- c(-0.30, -0.10, 0.10, 0.30)
FUNNEL_WIDTH <- 0.18
funnel_counts$x <- as.integer(funnel_counts$condition) +
  FUNNEL_DODGE[as.integer(funnel_counts$metric)]

# one total per bar, placed above the top of its stack
funnel_totals <- funnel_counts %>%
  group_by(condition, metric, x) %>%
  summarise(n = sum(n), .groups = "drop")

# In-segment labels for the ENUMERATE share only -- that is the quantity this
# revision exists to expose, and labelling all five values per condition would
# put a number on every mark. grep sits at the bottom, so the enumerate
# segment's midpoint is n_grep + n_enum/2. Dropped where the segment is too
# short to hold the text (CD7's 7 and 3 would overflow their own bar).
funnel_enum <- funnel_counts %>%
  filter(site_source %in% "enumerate") %>%
  left_join(funnel_counts %>%
              filter(site_source %in% "grep") %>%
              select(condition, metric, n_grep = n),
            by = c("condition", "metric")) %>%
  mutate(ymid = n_grep + n / 2)
funnel_enum <- funnel_enum[funnel_enum$n >= 0.06 * max(funnel_totals$n), ]

p_funnel <- ggplot(funnel_counts, aes(x = x, y = n, fill = fill)) +
  geom_col(position = position_stack(reverse = TRUE), width = FUNNEL_WIDTH,
           color = "#fcfcfb", linewidth = 0.5) +
  geom_text(data = funnel_totals, aes(x = x, y = n, label = format(n, big.mark = ",")),
            inherit.aes = FALSE, vjust = -0.45, size = 2.7) +
  geom_text(data = funnel_enum, aes(x = x, y = ymid, label = format(n, big.mark = ",")),
            inherit.aes = FALSE, size = 2.6, color = "#0b0b0b") +
  scale_fill_manual(values = funnel_palette, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = seq_along(conditions), labels = conditions,
                      expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                      labels = function(v) format(v, big.mark = ",")) +
  # 7 fills; 3 rows keeps the legend from eating the panel at this width
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  labs(x = NULL, y = "Peaks (bars 1-3) or edited genomic positions (bar 4)",
       title = "Differential peaks, GuideScan matches, edited R-loops and edited positions, by condition",
       subtitle = paste("Bars 1-3 are nested peak counts: every edited R-loop has a GuideScan match, and every match is a differential peak.",
                         "Bar 4 changes the unit to distinct edited genomic BASES, so it exceeds bar 3 wherever a peak carries >1 edit.",
                         "Stacked bars split by which finder located the site: grep (step 6, exact seed anchor; darker) vs enumerate",
                         "(step 7, no seed requirement; lighter). Within each bar the two segments are disjoint and sum to its total.",
                         sep = "\n")) +
  theme_minimal(base_size = 12, base_family = PLOT_FONT) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "peak_guidescan_edit_counts.pdf"), p_funnel, width = 13, height = 7)
ggsave(file.path(out_dir, "peak_guidescan_edit_counts.png"), p_funnel, width = 13, height = 7, dpi = 200)

for (i in seq_len(nrow(funnel_counts))) {
  cat(sprintf("  %-14s %-20s %-10s %6d (%.1f%% of peaks)\n",
              as.character(funnel_counts$condition[i]), funnel_counts$metric[i],
              ifelse(is.na(funnel_counts$site_source[i]), "-", funnel_counts$site_source[i]),
              funnel_counts$n[i], funnel_counts$pct_of_peaks[i]))
}

# =============================================================================
# SUBSTRATE-NORMALIZED editing rate.
#
# Reviewer objection this addresses: "peaks (R-loops) per edit" does not
# account for differences in substrate availability at off-target loci. If one
# condition's off-target protospacers simply contain fewer adenines in the
# editing window, it will show fewer edits per R-loop for reasons that have
# nothing to do with how efficiently the editor works -- and the comparison is
# outright invalid between editors with different substrates (ABE needs A, CBE
# needs C). The objection is correct, and it bites here: mean protospacer-strand
# adenines per off-target site is 3.62 / 3.60 / 3.37 / 4.67 for HEK4 /
# HEK4_noPalbo / HEK4_noSort / CD7 -- CD7 carries ~29% more substrate per site
# than HEK4_noSort, so part of CD7's favourable peaks-per-edit ratio is just
# having more adenines to edit.
#
# Fix: change the DENOMINATOR from "R-loops" to "available substrate bases".
# For an A>G edit at protospacer position p, that site's protospacer base at p
# must be an A, so the correct denominator at p is the number of off-target
# sites whose protospacer has an A at p:
#
#     rate(p) = A>G edits at p / sites with an A at p
#
# which is a true per-substrate probability in [0, 1] and is composition-free
# by construction. The same is done for T>C against available T. Summing
# numerator and denominator over p gives one substrate-normalized figure per
# condition (edits per 1,000 available adenines) that IS comparable across
# conditions -- and, unlike peaks-per-edit, would remain comparable across
# editors with different substrates.
#
# Off-target only (on-target excluded), matching the three mirrored bar plots.
#
# TWO KNOWN CAVEATS -- the denominator spans all 20 positions by choice, so
# both are live in the numbers below:
#
# 1. Positions 16-20 are SEED-INVARIANT. Stage A of 6A_guidescan.py requires an
#    exact --seed match, and with seed_offset 15 that pins spacer positions
#    16-20 to the seed's own bases in every single site: GGTGG for the HEK4
#    guide (0 adenines) but GGACC for CD7 (an A at position 18 in all 163
#    sites). So 163/762 = 21% of CD7's available_A is contributed by the search
#    method rather than by off-target biology, and CD7's pooled rate is
#    correspondingly diluted: 61.7 per 1k A over positions 1-20 vs 76.8 over
#    the variable positions 1-15. The three HEK4-locus conditions are
#    unaffected (their seed contributes no A). The per-position plot is immune
#    -- invariant positions simply have n_A = 0 and a rate of NA -- so this
#    only affects the pooled summary figure.
#
# 2. Composition is not independent of the guide. Off-target sites are found BY
#    HOMOLOGY to the spacer, so their base composition tracks the guide's:
#    the CD7 guide has 4 adenines vs HEK4's 2, which is most of why CD7 sites
#    average 4.67 A/site against HEK4's 3.62. Position 15 shows it directly --
#    78% A in HEK4 sites, 69% in CD7, where both guides carry an A, far above
#    the ~25% expected by chance.
#
# Consequence for interpretation: WITHIN a locus (HEK4 / HEK4_noPalbo /
# HEK4_noSort all share one guide and one seed) this normalization is clean and
# the conditions are directly comparable. ACROSS loci (any of those vs CD7)
# available substrate is set largely by guide composition and seed content, so
# treat the cross-locus ratio as indicative, not as a clean efficiency
# comparison.
# =============================================================================
substrate_by_position <- do.call(rbind, lapply(conditions, function(cond) {
  # BOTH site sources, matching the plots (see site_spacer_sequences: the
  # sequence is read from the reference over each site's own interval, so no
  # tool's match_sequence encoding is involved and step 6 and step 7 sites are
  # treated identically -- exactly as the plots treat them).
  sites <- read_sites_all(file.path(seed_dir, seed_files[[cond]]), cond)
  sites <- sites[sites$guidescan_match, ]
  on_t <- !is.na(sites$mismatches) & sites$mismatches == 0L &
    sites$rna_bulges == 0L & sites$dna_bulges == 0L
  sites <- sites[!on_t, ]
  spacer <- site_spacer_sequences(sites)
  spacer <- spacer[nchar(spacer) > 0]
  # Sites are NOT required to be exactly 20 bp. A bulged alignment's interval is
  # 20 - rna + dna, so it can be shorter or longer; dropping those would remove
  # them from the denominator while their edits stayed in the numerator and
  # would inflate the rate. Instead each site contributes only to the positions
  # it actually spans (short sites leave the tail NA, long ones are truncated at
  # 20), and colSums ignores the NAs.
  chars <- matrix(NA_character_, length(spacer), PROTOSPACER_LEN)
  for (i in seq_along(spacer)) {
    cs <- strsplit(spacer[i], "", fixed = TRUE)[[1]]
    k <- min(length(cs), PROTOSPACER_LEN)
    chars[i, seq_len(k)] <- cs[seq_len(k)]
  }
  data.frame(condition = cond, position = seq_len(PROTOSPACER_LEN),
             n_sites = length(spacer),
             n_A = as.integer(colSums(chars == "A", na.rm = TRUE)),
             n_T = as.integer(colSums(chars == "T", na.rm = TRUE)),
             stringsAsFactors = FALSE)
}))

# observed strand-resolved edits per position, off-target only
edits_by_position <- do.call(rbind, lapply(conditions, function(cond) {
  hits <- hits_by_condition[[cond]]
  hits <- hits[!hits$on_target, ]
  data.frame(condition = cond, position = seq_len(PROTOSPACER_LEN),
             # STRANDED: A>G only, whose substrate is a protospacer-strand A.
             # Aggregated: the pooled T>C class, whose substrate is an editable
             # adenine on EITHER strand -- hence the n_A + n_T denominator set
             # below rather than n_A alone.
             n_AG = if (STRANDED) as.integer(tabulate(hits$position[hits$ps_type == "A>G"], PROTOSPACER_LEN))
                     else as.integer(tabulate(hits$position[hits$snp_type == "T>C"], PROTOSPACER_LEN)),
             n_TC = if (STRANDED) as.integer(tabulate(hits$position[hits$ps_type == "T>C"], PROTOSPACER_LEN))
                     else 0L,
             stringsAsFactors = FALSE)
}))

substrate_rate <- merge(substrate_by_position, edits_by_position,
                         by = c("condition", "position"))
# Denominator: protospacer adenines when strand-resolved, A + T when pooled.
substrate_rate$n_denom <- if (STRANDED) substrate_rate$n_A else substrate_rate$n_A + substrate_rate$n_T
substrate_rate$rate_AG <- ifelse(substrate_rate$n_denom > 0, substrate_rate$n_AG / substrate_rate$n_denom, NA)
substrate_rate$rate_TC <- if (STRANDED) ifelse(substrate_rate$n_T > 0, substrate_rate$n_TC / substrate_rate$n_T, NA) else NA_real_
substrate_rate$condition <- factor(substrate_rate$condition, levels = conditions)
substrate_rate <- substrate_rate[order(substrate_rate$condition, substrate_rate$position), ]

write.csv(substrate_rate, file.path(out_dir, "substrate_normalized_rate.csv"), row.names = FALSE)

# ---- per-position A>G rate, one panel per condition ------------------------
rate_y_max <- max(substrate_rate$rate_AG, na.rm = TRUE)
rate_plots <- lapply(conditions, function(cond) {
  df <- substrate_rate[substrate_rate$condition == cond, ]
  ggplot(df, aes(x = position, y = rate_AG)) +
    geom_col(fill = snp_type_palette[["T>C"]], color = "black",
             linewidth = 0.2, width = 0.8, alpha = 0.9) +
    geom_text(aes(label = n_denom), vjust = -0.5, size = 2.2, color = "grey35") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    coord_cartesian(xlim = x_range, ylim = c(0, rate_y_max * 1.12)) +
    labs(x = X_LAB,
         y = if (STRANDED) "A>G edits per available adenine" else "T>C edits per editable adenine (either strand)",
         title = sprintf("%s (%d off-target sites)", cond, df$n_sites[1])) +
    theme_minimal(base_size = 11, base_family = PLOT_FONT) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", color = condition_palette[[cond]]))
})

p_rate <- wrap_plots(rate_plots, ncol = panel_ncol) +
  plot_annotation(
    title = "Substrate-normalized editing rate by protospacer position",
    subtitle = paste(
      "A>G edits (protospacer-strand adenine deamination) divided by the number of off-target sites that actually",
      "HAVE an adenine at that position -- a per-substrate probability, so it is not confounded by differences in",
      "adenine content between conditions. Small grey numbers are the denominators (available adenines). Shared y-axis.",
      sep = "\n"),
    theme = theme(text = element_text(family = PLOT_FONT), plot.title = element_text(face = "bold")))

ggsave(file.path(out_dir, "substrate_normalized_rate.pdf"), p_rate, width = panel_w, height = panel_h)
ggsave(file.path(out_dir, "substrate_normalized_rate.png"), p_rate, width = panel_w, height = panel_h, dpi = 200)

# ---- one normalized figure per condition, vs. the un-normalized ratio ------
substrate_summary <- substrate_rate %>%
  group_by(condition) %>%
  summarise(n_sites = first(n_sites),
            available_A = sum(n_A), available_T = sum(n_T),
            available_denom = sum(n_denom),
            edits_AG = sum(n_AG), edits_TC = sum(n_TC),
            .groups = "drop") %>%
  mutate(mean_denom_per_site = available_denom / n_sites,
         edits_per_1k_substrate = 1000 * edits_AG / available_denom,
         # the un-normalized metric, for side-by-side comparison
         edits_per_site = edits_AG / n_sites)
substrate_summary$condition <- factor(substrate_summary$condition, levels = conditions)

write.csv(substrate_summary, file.path(out_dir, "substrate_normalized_summary.csv"), row.names = FALSE)

p_norm <- ggplot(substrate_summary, aes(x = condition, y = edits_per_1k_substrate, fill = condition)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f\n(%s %s / %s edits)", edits_per_1k_substrate,
                                  format(available_denom, big.mark = ","),
                                  if (STRANDED) "A" else "A+T",
                                  format(edits_AG, big.mark = ","))),
            vjust = -0.25, size = 2.9, lineheight = 0.95) +
  scale_fill_manual(values = condition_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL,
       y = if (STRANDED) "A>G edits per 1,000 available adenines"
            else "T>C edits per 1,000 editable adenines (either strand)",
       title = paste0("Substrate-normalized editing efficiency, by condition",
                       if (STRANDED) "" else " (aggregated)"),
       subtitle = paste("Denominator is available substrate (protospacer-strand adenines across all off-target sites),",
                         "not R-loop count -- so this is not confounded by differing adenine content.", sep = "\n")) +
  theme_minimal(base_size = 12, base_family = PLOT_FONT) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "substrate_normalized_summary.pdf"), p_norm, width = 8, height = 6)
ggsave(file.path(out_dir, "substrate_normalized_summary.png"), p_norm, width = 8, height = 6, dpi = 200)

cat("\nSubstrate-normalized summary (off-target only):\n")
for (i in seq_len(nrow(substrate_summary))) {
  r <- substrate_summary[i, ]
  cat(sprintf("  %-14s sites=%5d  substrate=%6d (%.2f/site; A=%d T=%d)  edits=%4d  %6.1f per 1k %s   [unnormalized: %.3f edits/site]\n",
              as.character(r$condition), r$n_sites, r$available_denom, r$mean_denom_per_site,
              r$available_A, r$available_T, r$edits_AG, r$edits_per_1k_substrate,
              if (STRANDED) "A" else "A+T", r$edits_per_site))
}


# =============================================================================
# Protospacer edits as VCF  (folded in from Scripts/10_vcf.R)
#
# The analysis above answers "where in the spacer did editing happen" as tables
# and figures. This writes the same join back out as VCF, two files per
# condition per editor in <out-dir>/vcf/:
#
#   <cond>_observed_<TAG>_edits.vcf   signature edits called inside a scored
#                                     protospacer                 (NUMERATOR)
#   <cond>_possible_<TAG>_edits.vcf   every substrate base available to that
#                                     editor in a scored protospacer, edited
#                                     or not                    (DENOMINATOR)
#
# The second strictly contains the first (asserted, not assumed), so a pair is
# directly usable as a per-substrate rate. The CBE pass is not a second
# experiment -- every library here is ABE8e, so C>T at the same protospacers is
# the non-signature background, a within-sample specificity control.
#
# STRAND: an editor deaminates on the PROTOSPACER strand but VCF REF/ALT are
# reference-strand, and about half the scored sites are '-'. So ABE reads A>G on
# '+' and T>C on '-'. REF always matches the reference; INFO/PS_REF and PS_ALT
# carry the protospacer view, so the signature reads identically on every record.
#
# Written once, from the stranded pass only: the records are inherently
# strand-resolved (collapsing A>G with T>C is exactly what would destroy them),
# so the aggregated pass would emit the same bytes.
# =============================================================================
if (STRANDED) {

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

  # ==========================================================================
  # EDITOR-INDEPENDENT WORK -- done once per condition, shared by both passes.
  #
  # The sites, their reference sequences and the variant/site join do not
  # depend on which base the editor deaminates, so a second editor costs
  # essentially nothing beyond writing its two files.
  # ==========================================================================
  peaks <- read_sites_all(tsv_path, cond)
  sites <- peaks[peaks$guidescan_match, ]
  sites$on_target <- !is.na(sites$mismatches) & sites$mismatches == 0L &
    sites$rna_bulges == 0L & sites$dna_bulges == 0L
  sites$site_id <- seq_len(nrow(sites))
  if (sum(sites$on_target) != 1) {
    stop(sprintf("%s: expected exactly 1 on-target site, found %d",
                 cond, sum(sites$on_target)))
  }
  cat(sprintf("  %d peaks, %d with a scored protospacer (%d seed, %d enumerate), %d unscored\n",
              nrow(peaks), nrow(sites), sum(sites$site_source == "seed"),
              sum(sites$site_source == "enumerate"), sum(!peaks$guidescan_match)))

  site_label <- region_label(sites$site_chrom, sites$site_start, sites$site_end)
  peak_label <- region_label(sites$peak_chrom, sites$peak_start, sites$peak_end)

  # protospacer-strand sequence of every site, as a site x position character
  # matrix. NA past a short (RNA-bulged) site's end, truncated at 20 for a long
  # (DNA-bulged) one -- the same treatment substrate_by_position uses, which is
  # what makes the per-base counts equal 9_protospacer.R's available_* figures.
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
  # protospacer-oriented substitution: REF/ALT complemented on '-' sites. NOT
  # collapsed -- keeping A>G distinct from T>C (and C>T from G>A) is the entire
  # point, since the collapse is exactly what loses the edited strand.
  ps_ref <- ifelse(s_strand == "-", unname(COMPLEMENT[variants$ref[m]]), variants$ref[m])
  ps_alt <- ifelse(s_strand == "-", unname(COMPLEMENT[variants$alt[m]]), variants$alt[m])

  # ==========================================================================
  # PER-EDITOR PASSES
  # ==========================================================================
  for (ed in editors) {
    spec <- EDITOR_SPECS[[ed]]
    sig <- sprintf("%s>%s", spec$sub, spec$edit)
    cat(sprintf("\n  --- %s (protospacer-strand %s) ---\n", ed, sig))

    # ------------------------------------------------------------------
    # OBSERVED -- this editor's signature edits inside a scored protospacer
    # ------------------------------------------------------------------
    keep <- which(ps_ref == spec$sub & ps_alt == spec$edit)
    outside <- sum(ps_pos[keep] < 1L | ps_pos[keep] > PROTOSPACER_LEN)
    if (outside > 0) {
      warning(sprintf("%s/%s: %d %s edit(s) outside protospacer positions 1-%d (bulged site longer than the spacer); they cannot appear in the possible file",
                      cond, ed, outside, sig, PROTOSPACER_LEN))
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

    # ------------------------------------------------------------------
    # POSSIBLE -- every substrate base available to this editor
    # ------------------------------------------------------------------
    idx <- which(chars == spec$sub, arr.ind = TRUE)   # NA compares NA, so ignored
    idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
    p_site <- idx[, 1]
    p_pos  <- idx[, 2]

    plus <- sites$strand[p_site] == "+"
    # protospacer position p -> genomic coordinate, then the reference-strand
    # allele: a protospacer substrate base is itself on '+' and its complement
    # on '-' (ABE: A / T; CBE: C / G).
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

    # ------------------------------------------------------------------
    # CHECKS -- all fatal. 10_vcf.R also warned on drift from its own project's
    # published counts; those are dataset-specific and dropped here.
    # ------------------------------------------------------------------
    # 1. Every observed edit must be an available substrate base of the same
    #    site and position. This is the numerator-in-denominator property that
    #    makes the pair usable as a rate, and it cross-validates the strand
    #    handling: the two halves reach REF/ALT by different routes (Mutect2's
    #    call complemented, vs the reference sequence complemented), so a
    #    disagreement in either the complement or the position arithmetic
    #    breaks containment.
    stopifnot(all(okey %in% pkey))
    # 2. REF must be the actual hg38 base. The possible file is derived from the
    #    reference so this is near-tautological, but the observed file's REF comes
    #    from Mutect2, so this checks the join really landed on the site it claims.
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

# Side-by-side rate comparison, the reason BOTH is the default. On ABE8e
# libraries the CBE signature is not a second experiment -- it is the
# non-signature background at the very same protospacers, so the ratio of the
# two rates is a within-sample specificity estimate that needs no separate
# control library.
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

}
