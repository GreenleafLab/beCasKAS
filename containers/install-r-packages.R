# Install Bioconductor packages for differential_peaks step

install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)

bioc_pkgs <- c("DiffBind", "DESeq2", "edgeR", "GenomicRanges", "rtracklayer")
cran_pkgs <- c("optparse", "ggplot2", "dplyr", "tidyr", "scales")

BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

# Verify
loaded <- sapply(c(bioc_pkgs, cran_pkgs), function(p) {
    requireNamespace(p, quietly = TRUE)
})
if (!all(loaded)) {
    stop("Some packages failed to install: ", paste(names(loaded)[!loaded], collapse = ", "))
}
cat("All R packages installed successfully.\n")
