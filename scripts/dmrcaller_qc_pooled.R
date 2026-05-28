#!/usr/bin/env Rscript
.libPaths("~/R/library")
# dmrcaller_qc_pooled.R
# Coverage and Correlation QC — all 12 samples pooled into one dataset
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/qc/for_meeting"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- c(
  "Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
  "Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3",
  "ASO_CTRL_1",      "ASO_CTRL_2",      "ASO_CTRL_3",
  "ASO_VPA_1",       "ASO_VPA_2",       "ASO_VPA_3"
)

CHROM     <- "chr1"
breaks    <- c(1, 5, 10, 15, 20, 30, 50)
distances <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)

# Load all 12 samples
message("Loading all 12 samples...")
meth_list <- list()
for (sample in SAMPLES) {
  path <- file.path(COV_DIR, paste0(sample, "_", CHROM, ".CpG_report.txt.gz"))
  if (!file.exists(path)) { message("Missing: ", path); next }
  message("  Loading: ", sample)
  meth_list[[sample]] <- readBismark(path)
}

# Pool all into one dataset
message("Pooling all 12 samples...")
methylationDataAll <- poolMethylationDatasets(GRangesList(meth_list))
message("Pooled.")

# Coverage
message("Plotting coverage...")
pdf(file.path(OUT_DIR, "final_coverage_pooled.pdf"), width = 9, height = 6)
par(mar = c(4, 4, 3, 1) + 0.1)
cov <- computeMethylationDataCoverage(methylationDataAll, context = "CG", breaks = breaks)
plot(breaks, cov, type = "l", col = "#1D6FA4", lwd = 2.5,
     ylim = c(0, 1), xlab = "Coverage depth", ylab = "Proportion of CpGs",
     main = "CpG Coverage — All 12 Samples Pooled (chr1)")
abline(v = c(5, 10), lty = 2, col = "grey50")
legend("topright", legend = "All samples pooled (n=12)",
       col = "#1D6FA4", lwd = 2.5, bty = "n")
dev.off()
message("Saved: final_coverage_pooled.pdf")

# Spatial correlation
message("Plotting spatial correlation...")
pdf(file.path(OUT_DIR, "final_spatial_correlation_pooled.pdf"), width = 9, height = 6)
par(mar = c(4, 4, 3, 1) + 0.1)
cor <- computeMethylationDataSpatialCorrelation(
  methylationDataAll, context = "CG", distances = distances)
plot(distances, cor, type = "l", col = "#1D6FA4", lwd = 2.5,
     ylim = c(0, 1), xlab = "Distance (bp)", ylab = "Correlation",
     main = "CpG Spatial Correlation — All 12 Samples Pooled (chr1)",
     log = "x")
dev.off()
message("Saved: final_spatial_correlation_pooled.pdf")

message("All done. Outputs in: ", OUT_DIR)
