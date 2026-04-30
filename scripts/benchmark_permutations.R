.libPaths("~/R/library")
# benchmark_permutations.R
# Multiple permutation benchmark — runs stratified scramble 20 times
# with different random seeds to get mean +/- SD signal/noise
# Tests bins strict at all window sizes on chr1
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 248956422
N_PERMS    <- 20
SEEDS      <- 1:N_PERMS

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Loading chr1 data...")
aso_vpa <- readBismarkPool(c(
  file.path(COV_DIR, "ASO_VPA_1_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_2_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_3_chr1.CpG_report.txt.gz")
))
aso_ctrl <- readBismarkPool(c(
  file.path(COV_DIR, "ASO_CTRL_1_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_CTRL_2_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_CTRL_3_chr1.CpG_report.txt.gz")
))
message("ASO_VPA CpGs:  ", length(aso_vpa))
message("ASO_CTRL CpGs: ", length(aso_ctrl))

scramble_stratified <- function(dat1, dat2, seed=42) {
  set.seed(seed)
  coverage <- mcols(dat1)$readsN
  strata <- cut(coverage,
                breaks = c(0, 5, 10, 20, 50, Inf),
                labels = c("1-5", "6-10", "11-20", "21-50", "50+"),
                include.lowest = TRUE)
  idx <- seq_along(coverage)
  for (s in levels(strata)) {
    stratum_idx <- which(strata == s)
    if (length(stratum_idx) > 1) idx[stratum_idx] <- sample(stratum_idx)
  }
  dat1_scr <- dat1
  dat2_scr <- dat2
  mcols(dat1_scr)$readsM <- mcols(dat1)$readsM[idx]
  mcols(dat1_scr)$readsN <- mcols(dat1)$readsN[idx]
  mcols(dat2_scr)$readsM <- mcols(dat2)$readsM[idx]
  mcols(dat2_scr)$readsN <- mcols(dat2)$readsN[idx]
  list(dat1=dat1_scr, dat2=dat2_scr)
}

window_sizes <- c(100, 200, 300, 500, 1000, 2000)
results <- list()

message("\nComputing real DMR counts...")
for (ws in window_sizes) {
  regions <- GRanges("chr1", IRanges(
    start = seq(1, REGION_END, by = ws), width = ws))
  dmrs_real <- computeDMRs(aso_vpa, aso_ctrl,
    regions=regions, context="CG", method="bins", binSize=ws,
    test="fisher", pValueThreshold=0.01,
    minCytosinesCount=4, minProportionDifference=0.2,
    minGap=0, minSize=50, minReadsPerCytosine=4,
    cores=32, parallel=TRUE)
  n_real <- length(dmrs_real)
  message("  ws=", ws, " real DMRs: ", n_real)
  results[[paste("real", ws)]] <- data.frame(
    window_size=ws, seed=NA, type="real", n_dmrs=n_real)
}

for (seed in SEEDS) {
  message("\n--- Permutation ", seed, "/", N_PERMS, " ---")
  scr <- scramble_stratified(aso_vpa, aso_ctrl, seed=seed)
  for (ws in window_sizes) {
    regions <- GRanges("chr1", IRanges(
      start = seq(1, REGION_END, by = ws), width = ws))
    dmrs_scr <- tryCatch({
      computeDMRs(scr$dat1, scr$dat2,
        regions=regions, context="CG", method="bins", binSize=ws,
        test="fisher", pValueThreshold=0.01,
        minCytosinesCount=4, minProportionDifference=0.2,
        minGap=0, minSize=50, minReadsPerCytosine=4,
        cores=32, parallel=TRUE)
    }, error=function(e) NULL)
    n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
    message("  ws=", ws, " seed=", seed, " scrambled: ", n_scr)
    results[[paste("scr", ws, seed)]] <- data.frame(
      window_size=ws, seed=seed, type="scrambled", n_dmrs=n_scr)
  }
}

df <- do.call(rbind, results)
write.csv(df, file.path(OUT_DIR, "benchmark_permutations_raw.csv"), row.names=FALSE)

real_counts <- df[df$type == "real", c("window_size", "n_dmrs")]
names(real_counts)[2] <- "n_real"

scr_sub <- df[df$type == "scrambled",]
scr_mean <- aggregate(n_dmrs ~ window_size, data=scr_sub, FUN=mean, na.rm=TRUE)
scr_sd   <- aggregate(n_dmrs ~ window_size, data=scr_sub, FUN=sd, na.rm=TRUE)
names(scr_mean)[2] <- "mean_scrambled"
names(scr_sd)[2]   <- "sd_scrambled"

summary_df <- merge(real_counts, scr_mean, by="window_size")
summary_df <- merge(summary_df, scr_sd, by="window_size")
summary_df$ratio_mean  <- round(summary_df$n_real / summary_df$mean_scrambled, 3)
summary_df$ratio_lower <- round(summary_df$n_real / (summary_df$mean_scrambled + summary_df$sd_scrambled), 3)
summary_df$ratio_upper <- round(summary_df$n_real / (summary_df$mean_scrambled - summary_df$sd_scrambled), 3)

write.csv(summary_df,
          file.path(OUT_DIR, "benchmark_permutations_summary.csv"),
          row.names=FALSE)

message("\n=== PERMUTATION RESULTS ===")
print(summary_df, row.names=FALSE)
message("\nSaved to: ", OUT_DIR)
