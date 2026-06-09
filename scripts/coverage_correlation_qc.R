#!/usr/bin/env Rscript
.libPaths("~/R/library")
library(DMRcaller)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/qc/for_meeting"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- c("Scramble_CTRL", "Scramble_VPA", "ASO_CTRL", "ASO_VPA")

reps <- 1:3

CHROM     <- "chr1"
breaks    <- c(1, 5, 10, 15, 20, 30, 50)
distances <- c(5, 10, 20, 50, 100, 200, 500)

message("Loading all 12 samples...")
meth_list <- list()
for (sample in SAMPLES) {
  meth_list[[sample]] <- list()
  for(rep in reps){
    path <- file.path(COV_DIR, paste0(sample, "_", rep,"_", CHROM, ".CpG_report.txt.gz"))
    if (!file.exists(path)) { message("Missing: ", path); next }
    message("  Loading: ", sample)
    meth_list[[sample]][[rep]] <- readBismark(path)
  }
}

message("Pooling all 12 samples...")
methylationDataPooled <- list()
coverage <- list()
for (sample in SAMPLES) {
  methylationDataPooled[[sample]] <- poolMethylationDatasets(GRangesList(meth_list[[sample]]))
  message(paste0("Pooled. ", sample))
  coverage[[sample]] <- computeMethylationDataCoverage(methylationDataPooled[[sample]], context = "CG", breaks = breaks)
}

message("Pooling all 4 conditions into combined dataset...")
methylationDataAll <- poolMethylationDatasets(GRangesList(methylationDataPooled))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

message("Plotting coverage...")
pdf(file.path(OUT_DIR, "radu_coverage_pooled.pdf"), width = 9, height = 6)
par(mar = c(4, 4, 3, 1) + 0.1)
cov <- computeMethylationDataCoverage(methylationDataAll, context = "CG", breaks = breaks)
plot(breaks, coverage[[1]], type = "n", col = "#1D6FA4", lwd = 2.5,
     ylim = c(0, 1), xlab = "Coverage depth", ylab = "Proportion of CpGs",
     main = "CpG Coverage — All 12 Samples Pooled (chr1)")
for (i in 1:length(SAMPLES)) {
  lines(breaks, coverage[[i]], type="l", col=cbbPalette[i], lty=1, lwd=2)
}
abline(v = c(5, 10), lty = 2, col = "grey50")
legend("topright", legend = SAMPLES,
       col = cbbPalette[1:length(SAMPLES)], lwd = 2.5, bty = "n")
dev.off()
message("Saved: radu_coverage_pooled.pdf")

message("Plotting spatial correlation...")
pdf(file.path(OUT_DIR, "radu_spatial_correlation_pooled.pdf"), width = 9, height = 6)
par(mar = c(4, 4, 3, 1) + 0.1)
cor <- computeMethylationDataSpatialCorrelation(
  methylationDataAll, context = "CG", distances = distances)
plot(distances, cor, type = "l", col = "#1D6FA4", lwd = 2.5,
     ylim = c(0, 1), xlab = "Distance (bp)", ylab = "Correlation",
     main = "CpG Spatial Correlation — All 12 Samples Pooled (chr1)",
     log = "x")
dev.off()
message("Saved: radu_spatial_correlation_pooled.pdf")

message("All done. Outputs in: ", OUT_DIR)
