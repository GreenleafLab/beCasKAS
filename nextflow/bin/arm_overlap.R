#!/usr/bin/env Rscript
# Overlap between the arms' variant-calling intervals for one case: a table over every
# combination of arms, and a Venn diagram of the same counts.
#
# usage: arm_overlap.R <case_id> <arm>=<bed> <arm>=<bed> [<arm>=<bed>]
# An empty BED (or /dev/null) is an arm that ran and found nothing. Writes arm_overlap.tsv and
# arm_overlap_venn.pdf to the working directory.
#
#   regions: merged intervals (all arms pooled, merged like bedtools merge) whose contributing
#            arms are exactly this combination
#   bp:      bases covered by exactly this combination of arms

suppressPackageStartupMessages(library(GenomicRanges))
options(scipen = 99)

args    <- commandArgs(trailingOnly = TRUE)
case_id <- args[1]
arms    <- sub("=.*", "", args[-1])
beds    <- setNames(sub("^[^=]*=", "", args[-1]), arms)
if (!length(arms) %in% 2:3) stop("arm_overlap.R draws 2 or 3 arms, got ", length(arms), call. = FALSE)

read_bed <- function(f) {
  if (file.size(f) == 0) return(GRanges())
  d <- read.table(f, sep = "\t", quote = "", comment.char = "#", colClasses = "character")
  reduce(GRanges(d[[1]], IRanges(as.integer(d[[2]]) + 1L, as.integer(d[[3]]))))
}
sets <- lapply(beds, read_bed)

# "101": overlaps arms 1 and 3, not 2
membership <- function(gr) {
  m <- do.call(cbind, lapply(sets, function(s) overlapsAny(gr, s)))
  if (length(gr)) apply(m * 1L, 1, paste, collapse = "") else character()
}
pooled   <- suppressWarnings(do.call(c, unname(sets)))
clusters <- membership(reduce(pooled))
pieces   <- disjoin(pooled)
piece_of <- membership(pieces)

combos <- expand.grid(rep(list(1:0), length(arms)))
names(combos) <- arms
combos <- combos[rowSums(combos) > 0, , drop = FALSE]
keys   <- apply(combos, 1, paste, collapse = "")
o      <- order(rowSums(combos), -strtoi(keys, base = 2))   # singles first, in arm order
combos <- combos[o, , drop = FALSE]
keys   <- keys[o]

tab <- data.frame(case = case_id, combos,
                  regions = vapply(keys, function(k) sum(clusters == k), integer(1)),
                  bp      = vapply(keys, function(k) sum(as.numeric(width(pieces)[piece_of == k])), numeric(1)),
                  row.names = NULL, check.names = FALSE)
write.table(tab, "arm_overlap.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(tab, row.names = FALSE)

fmt_bp <- function(b) if (b >= 1e6) sprintf("%.1f Mb", b / 1e6) else if (b >= 1e3) sprintf("%.1f kb", b / 1e3) else sprintf("%d bp", as.integer(b))

# Fixed layouts: circle centres, arm-name positions, and one label position per combination key
if (length(arms) == 2) {
  ctr  <- rbind(c(-0.6, 0), c(0.6, 0))
  name <- rbind(c(-0.6, 1.25), c(0.6, 1.25))
  at   <- list("10" = c(-1, 0), "01" = c(1, 0), "11" = c(0, 0))
} else {
  ctr  <- rbind(c(-0.6, 0.4), c(0.6, 0.4), c(0, -0.6))
  name <- rbind(c(-0.9, 1.65), c(0.9, 1.65), c(0, -1.85))
  at   <- list("100" = c(-1.05, 0.75), "010" = c(1.05, 0.75), "001" = c(0, -1.15),
               "110" = c(0, 0.95), "101" = c(-0.7, -0.4), "011" = c(0.7, -0.4), "111" = c(0, 0.05))
}
cols <- c("#1565C0", "#D97827", "#0B6623")[seq_along(arms)]

pdf("arm_overlap_venn.pdf", width = 6, height = 6)
par(mar = c(1, 1, 3, 1))
plot.new()
plot.window(c(-2.1, 2.1), c(-2.1, 2.1), asp = 1)
symbols(ctr[, 1], ctr[, 2], circles = rep(1, length(arms)), inches = FALSE, add = TRUE,
        fg = cols, bg = adjustcolor(cols, 0.15), lwd = 2)
text(name[, 1], name[, 2], arms, col = cols, font = 2)
for (i in seq_along(keys)) {
  xy <- at[[keys[i]]]
  text(xy[1], xy[2], sprintf("%d\n%s", tab$regions[i], fmt_bp(tab$bp[i])), cex = 0.8)
}
title(case_id)
mtext("merged regions / bp covered by exactly those arms", side = 1, cex = 0.8)
invisible(dev.off())
