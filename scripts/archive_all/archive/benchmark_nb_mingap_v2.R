.libPaths("~/R/library")
library(DMRcaller)
COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 248956422
N_SEEDS    <- 20
CHECKPOINT_FILE <- file.path(OUT_DIR, "checkpoint_nb_mingap_v2.rds")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

mingap_values  <- c(100, 200, 300, 500, 1000, 2000)
minsize_values <- c(100, 200)

scramble_data <- function(dat1, dat2, seed=42) {
  set.seed(seed)
  idx1 <- sample(seq_len(length(dat1)), replace=FALSE)
  idx2 <- sample(seq_len(length(dat2)), replace=FALSE)
  dat1_scr <- dat1; dat2_scr <- dat2
  mcols(dat1_scr)$readsM <- mcols(dat1)$readsM[idx1]
  mcols(dat1_scr)$readsN <- mcols(dat1)$readsN[idx1]
  mcols(dat2_scr)$readsM <- mcols(dat2)$readsM[idx2]
  mcols(dat2_scr)$readsN <- mcols(dat2)$readsN[idx2]
  list(dat1=dat1_scr, dat2=dat2_scr)
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

region <- GRanges(CHROM, IRanges(1, REGION_END))

results <- list()
completed_keys <- character(0)
if (file.exists(CHECKPOINT_FILE)) {
  message("Loading checkpoint: ", CHECKPOINT_FILE)
  ck <- readRDS(CHECKPOINT_FILE)
  results        <- ck$results
  completed_keys <- ck$completed_keys
  message("  Completed: ", paste(completed_keys, collapse=", "))
} else {
  message("Starting from scratch.")
}

save_checkpoint <- function() {
  saveRDS(list(results=results, completed_keys=completed_keys), CHECKPOINT_FILE)
  message("  Checkpoint saved.")
}

for (mg in mingap_values) {
  for (ms in minsize_values) {
    key <- paste(mg, ms, sep="_")
    if (key %in% completed_keys) {
      message("\n--- Skipping minGap=", mg, " minSize=", ms, " ---"); next
    }
    message("\n--- minGap=", mg, " | minSize=", ms, " | Real ---")
    dmrs_real <- computeDMRs(aso_vpa, aso_ctrl,
      regions=region, context="CG", method="neighbourhood",
      test="score", pValueThreshold=0.01,
      minCytosinesCount=4, minProportionDifference=0.2,
      minGap=mg, minSize=ms,
      minReadsPerCytosine=4, cores=16, parallel=TRUE)
    n_real <- length(dmrs_real)
    message("  Real DMRs: ", n_real)

    message("--- ", N_SEEDS, " permutations ---")
    scr_counts <- integer(N_SEEDS)
    for (s in seq_len(N_SEEDS)) {
      scr_s <- scramble_data(aso_vpa, aso_ctrl, seed=s)
      d_s <- computeDMRs(scr_s$dat1, scr_s$dat2,
        regions=region, context="CG", method="neighbourhood",
        test="score", pValueThreshold=0.01,
        minCytosinesCount=4, minProportionDifference=0.2,
        minGap=mg, minSize=ms,
        minReadsPerCytosine=4, cores=16, parallel=TRUE)
      scr_counts[s] <- length(d_s)
      message("    seed ", s, "/", N_SEEDS, ": ", scr_counts[s])
    }
    mean_scr <- mean(scr_counts)
    sd_scr   <- sd(scr_counts)
    ratio    <- ifelse(mean_scr==0, Inf, round(n_real/mean_scr, 3))
    message("  mean +/- sd: ", round(mean_scr,1), " +/- ", round(sd_scr,1))
    message("  S/N: ", ratio)

    results[[key]] <- data.frame(
      method="neighbourhood_mingap", mingap=mg, minsize=ms,
      scramble_method="independent_20seeds",
      n_real=n_real, mean_scrambled=mean_scr,
      sd_scrambled=sd_scr, ratio=ratio)
    completed_keys <- c(completed_keys, key)
    save_checkpoint()
  }
}

df <- do.call(rbind, results)
out_csv <- file.path(OUT_DIR, "benchmark_nb_mingap_v2.csv")
write.csv(df, out_csv, row.names=FALSE)
message("\n=== RESULTS ===")
print(df, row.names=FALSE)
message("Saved: ", out_csv)
