.libPaths("~/R/library")
# benchmark_permutations.R
# Multiple permutation benchmark with resume from checkpoint
# Random (non-stratified) shuffling as null model
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 248956422
N_PERMS    <- 20
SEEDS      <- 1:N_PERMS
CHECKPOINT_FILE <- file.path(OUT_DIR, "benchmark_permutations_checkpoint.rds")

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

# -------------------------------------------------------------------
# Radu's fully random (non-stratified) permutation
# -------------------------------------------------------------------
scramble_data <- function(dat1, dat2, seed = 42) {
  set.seed(seed)
  idx <- sample(seq_len(length(dat1)), size = length(dat1), replace = FALSE)
  dat1_scr <- dat1
  dat2_scr <- dat2
  mcols(dat1_scr)$readsM <- mcols(dat1)$readsM[idx]
  mcols(dat1_scr)$readsN <- mcols(dat1)$readsN[idx]
  mcols(dat2_scr)$readsM <- mcols(dat2)$readsM[idx]
  mcols(dat2_scr)$readsN <- mcols(dat2)$readsN[idx]
  list(dat1 = dat1_scr, dat2 = dat2_scr)
}

window_sizes <- c(100, 200, 300, 500, 1000, 2000)

# -------------------------------------------------------------------
# Load or initialize checkpoint
# -------------------------------------------------------------------
results <- list()
completed_real <- FALSE
completed_seeds <- integer(0)

if (file.exists(CHECKPOINT_FILE)) {
  message("Loading checkpoint from ", CHECKPOINT_FILE)
  ck <- readRDS(CHECKPOINT_FILE)
  results         <- ck$results
  completed_real  <- isTRUE(ck$completed_real)
  completed_seeds <- ck$completed_seeds
  message("  completed_real: ", completed_real)
  if (length(completed_seeds) > 0) {
    message("  completed_seeds: ", paste(completed_seeds, collapse = ", "))
  }
} else {
  message("No checkpoint found; starting from scratch.")
}

save_checkpoint <- function() {
  ck <- list(
    results         = results,
    completed_real  = completed_real,
    completed_seeds = completed_seeds
  )
  saveRDS(ck, CHECKPOINT_FILE)
  message("Checkpoint saved to ", CHECKPOINT_FILE)
}

# -------------------------------------------------------------------
# Real DMR counts (once per window size)
# -------------------------------------------------------------------
if (!completed_real) {
  message("\nComputing real DMR counts...")
  for (ws in window_sizes) {
    key <- paste("real", ws, sep = "_")
    if (!is.null(results[[key]])) {
      next
    }
    regions <- GRanges(CHROM, IRanges(
      start = seq(1, REGION_END, by = ws), width = ws))
    dmrs_real <- computeDMRs(aso_vpa, aso_ctrl,
      regions=regions, context="CG", method="bins", binSize=ws,
      test="fisher", pValueThreshold=0.01,
      minCytosinesCount=4, minProportionDifference=0.2,
      minGap=0, minSize=50, minReadsPerCytosine=4,
      cores=32, parallel=TRUE)
    n_real <- length(dmrs_real)
    message("  ws=", ws, " real DMRs: ", n_real)
    results[[key]] <- data.frame(
      window_size=ws, seed=NA, type="real", n_dmrs=n_real)
  }
  completed_real <- TRUE
  save_checkpoint()
} else {
  message("Real DMR counts already completed; skipping.")
}

# -------------------------------------------------------------------
# Permutation loop with resume
# -------------------------------------------------------------------
for (seed in SEEDS) {
  if (seed %in% completed_seeds) {
    message("\n--- Seed ", seed, " already completed; skipping ---")
    next
  }
  message("\n--- Permutation ", seed, "/", N_PERMS, " ---")
  scr <- tryCatch({
    scramble_data(aso_vpa, aso_ctrl, seed=seed)
  }, error = function(e) {
    message("  Error in scramble_data for seed ", seed, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(scr)) {
    next
  }
  for (ws in window_sizes) {
    key <- paste("scr", ws, seed, sep = "_")
    if (!is.null(results[[key]])) {
      next
    }
    regions <- GRanges(CHROM, IRanges(
      start = seq(1, REGION_END, by = ws), width = ws))
    dmrs_scr <- tryCatch({
      computeDMRs(scr$dat1, scr$dat2,
        regions=regions, context="CG", method="bins", binSize=ws,
        test="fisher", pValueThreshold=0.01,
        minCytosinesCount=4, minProportionDifference=0.2,
        minGap=0, minSize=50, minReadsPerCytosine=4,
        cores=32, parallel=TRUE)
    }, error=function(e) {
      message("  Error in computeDMRs, ws=", ws, " seed=", seed, ": ", conditionMessage(e))
      NULL
    })
    n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
    message("  ws=", ws, " seed=", seed, " scrambled: ", n_scr)
    results[[key]] <- data.frame(
      window_size=ws, seed=seed, type="scrambled", n_dmrs=n_scr)
  }
  completed_seeds <- sort(unique(c(completed_seeds, seed)))
  save_checkpoint()
}

# -------------------------------------------------------------------
# Combine and write outputs
# -------------------------------------------------------------------
message("\nCombining results...")
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
