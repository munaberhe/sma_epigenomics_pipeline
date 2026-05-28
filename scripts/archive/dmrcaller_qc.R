
.libPaths("~/R/library")
# dmrcaller_qc.R
# Coverage and Correlation QC using DMRcaller native functions
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/qc/dmrcaller"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- list(
  Scramble_CTRL = c("Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3"),
  Scramble_VPA  = c("Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3"),
  ASO_CTRL      = c("ASO_CTRL_1",      "ASO_CTRL_2",      "ASO_CTRL_3"),
  ASO_VPA       = c("ASO_VPA_1",       "ASO_VPA_2",       "ASO_VPA_3")
)

CHROM <- "chr1"
breaks <- c(1, 5, 10, 15, 20, 30, 50)
distances <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)

group_colours <- c(
  Scramble_CTRL = "#065A82",
  Scramble_VPA  = "#1C7293",
  ASO_CTRL      = "#02C39A",
  ASO_VPA       = "#F59E0B"
)

message("Loading chr1 data for all samples...")
meth_data <- list()
for (group in names(SAMPLES)) {
  for (sample in SAMPLES[[group]]) {
    path <- file.path(COV_DIR, paste0(sample, "_", CHROM, ".CpG_report.txt.gz"))
    if (!file.exists(path)) { message("Missing: ", path); next }
    message("  Loading: ", sample)
    dat <- readBismark(path)
    dat <- dat[dat$readsN >= 1]
    meth_data[[sample]] <- dat
  }
}
message("Loaded ", length(meth_data), " samples")

# Coverage per group
message("Computing coverage...")
pdf(file.path(OUT_DIR, "dmrcaller_coverage.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (group in names(SAMPLES)) {
  cov_list <- list()
  for (sample in SAMPLES[[group]]) {
    if (!sample %in% names(meth_data)) next
    cov_list[[sample]] <- computeMethylationDataCoverage(
      meth_data[[sample]], context = "CG", breaks = breaks)
  }
  if (length(cov_list) == 0) next
  cov_mat <- do.call(rbind, cov_list)
  plot(breaks, cov_mat[1,], type="l", col=group_colours[group],
       ylim=c(0,1), xlab="Coverage depth", ylab="Proportion of CpGs",
       main=paste0("Coverage — ", group), lwd=2)
  for (i in seq_len(nrow(cov_mat)))
    lines(breaks, cov_mat[i,], col=group_colours[group], lwd=1.5, lty=i)
  abline(v=c(5,10), lty=2, col="grey50")
  legend("topright", legend=SAMPLES[[group]], lty=1:3,
         col=rep(group_colours[group],3), cex=0.7)
}
dev.off()
message("Saved: dmrcaller_coverage.pdf")

# Cross-group coverage comparison
pdf(file.path(OUT_DIR, "dmrcaller_coverage_comparison.pdf"), width = 10, height = 7)
par(mar = c(4, 4, 3, 1))
plot(NULL, xlim=c(1,50), ylim=c(0,1),
     xlab="Coverage depth", ylab="Proportion of CpGs",
     main="Coverage Comparison — All Groups (chr1)")
abline(v=c(5,10), lty=2, col="grey80")
lty_i <- 1
for (group in names(SAMPLES)) {
  for (sample in SAMPLES[[group]]) {
    if (!sample %in% names(meth_data)) next
    cov <- computeMethylationDataCoverage(
      meth_data[[sample]], context="CG", breaks=breaks)
    lines(breaks, cov, col=group_colours[group], lwd=1.5, lty=lty_i)
  }
  lty_i <- lty_i + 1
}
legend("bottomleft", legend=names(group_colours), col=group_colours, lwd=2, cex=0.85)
dev.off()
message("Saved: dmrcaller_coverage_comparison.pdf")

# Spatial correlation per group
message("Computing spatial correlation...")
pdf(file.path(OUT_DIR, "dmrcaller_spatial_correlation.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (group in names(SAMPLES)) {
  cor_list <- list()
  for (sample in SAMPLES[[group]]) {
    if (!sample %in% names(meth_data)) next
    message("  Spatial correlation: ", sample)
    cor_list[[sample]] <- computeMethylationDataSpatialCorrelation(
      meth_data[[sample]], context = "CG", distances = distances)
  }
  if (length(cor_list) == 0) next
  plot(distances, cor_list[[1]], type="l", col=group_colours[group],
       ylim=c(0,1), xlab="Distance (bp)", ylab="Correlation",
       main=paste0("Spatial Correlation — ", group), log="x", lwd=2)
  for (i in seq_along(cor_list))
    lines(distances, cor_list[[i]], col=group_colours[group], lwd=1.5, lty=i)
  legend("topright", legend=SAMPLES[[group]], lty=1:3,
         col=rep(group_colours[group],3), cex=0.7)
}
dev.off()
message("Saved: dmrcaller_spatial_correlation.pdf")

message("\nDone. QC outputs in: ", OUT_DIR)
