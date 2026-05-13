#!/usr/bin/env Rscript
.libPaths("~/R/library")
# benchmark_nb_mingap_radu.R
# Neighbourhood minGap sweep using Radu's random permutation null model
# With checkpoint saving after each run
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 248956422
CHECKPOINT_FILE <- file.path(OUT_DIR, "checkpoint_nb_mingap_radu.rds")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

mingap_values <- c(100, 200, 300, 500, 1000, 2000)
minsize_values <- c(100, 200)

# Radu's fully random permutation
scramble_data <- function(dat1, dat2, seed = 42) {
  set.seed(seed)
  idx <- sample(seq_len(length(dat1)), size = length(dat1), replace = FALSE)
  dat1_scr <- dat1; dat2_scr <- dat2
  mcols(dat1_scr)$readsM <- mcols(dat1)$readsM[idx]
  mcols(dat1_scr)$readsN <- mcols(dat1)$readsN[idx]
  mcols(dat2_scr)$readsM <- mcols(dat2)$readsM[idx]
  mcols(dat2_scr)$readsN <- mcols(dat2)$readsN[idx]
  list(dat1 = dat1_scr, dat2 = dat2_scr)
}

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
message("Data loaded.")

message("Scrambling data (seed=42)...")
scr <- scramble_data(aso_vpa, aso_ctrl, seed = 42)

region <- GRanges(CHROM, IRanges(1, REGION_END))

# Load checkpoint if exists
results <- list()
completed_keys <- character(0)

if (file.exists(CHECKPOINT_FILE)) {
  message("Loading checkpoint from ", CHECKPOINT_FILE)
  ck <- readRDS(CHECKPOINT_FILE)
  results        <- ck$results
  completed_keys <- ck$completed_keys
  message("  Completed runs: ", paste(completed_keys, collapse=", "))
} else {
  message("No checkpoint found; starting from scratch.")
}

save_checkpoint <- function() {
  saveRDS(list(results=results, completed_keys=completed_keys), CHECKPOINT_FILE)
  message("  Checkpoint saved.")
}

for (mg in mingap_values) {
  for (ms in minsize_values) {
    key <- paste(mg, ms, sep="_")
    if (key %in% completed_keys) {
      message("\n--- Skipping minGap=", mg, " | minSize=", ms, " (already done) ---")
      next
    }

    message("\n--- radu | minGap=", mg, " | minSize=", ms, " | Real ---")
    dmrs_real <- computeDMRs(aso_vpa, aso_ctrl,
      regions=region, context="CG", method="neighbourhood",
      test="score", pValueThreshold=0.01,
      minCytosinesCount=4, minProportionDifference=0.2,
      minGap=mg, minSize=ms,
      minReadsPerCytosine=4, cores=16, parallel=TRUE)
    n_real <- length(dmrs_real)
    message("  Real DMRs: ", n_real)

    message("--- radu | minGap=", mg, " | minSize=", ms, " | Scrambled ---")
    dmrs_scr <- computeDMRs(scr$dat1, scr$dat2,
      regions=region, context="CG", method="neighbourhood",
      test="score", pValueThreshold=0.01,
      minCytosinesCount=4, minProportionDifference=0.2,
      minGap=mg, minSize=ms,
      minReadsPerCytosine=4, cores=16, parallel=TRUE)
    n_scr <- length(dmrs_scr)
    ratio <- ifelse(n_scr == 0, Inf, round(n_real / n_scr, 3))
    message("  Scrambled DMRs: ", n_scr)
    message("  S/N: ", ratio)

    results[[key]] <- data.frame(
      method="neighbourhood_mingap", mingap=mg, minsize=ms,
      scramble_method="radu",
      n_real=n_real, n_scrambled=n_scr, ratio=ratio)

    completed_keys <- c(completed_keys, key)
    save_checkpoint()
  }
}

df <- do.call(rbind, results)
out_csv <- file.path(OUT_DIR, "benchmark_nb_mingap_radu.csv")
write.csv(df, out_csv, row.names=FALSE)

message("\n=== NEIGHBOURHOOD minGap RESULTS (Radu permutation) ===")
print(df, row.names=FALSE)
message("\nSaved to: ", out_csv)
