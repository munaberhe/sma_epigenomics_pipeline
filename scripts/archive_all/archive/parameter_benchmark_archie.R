.libPaths("~/R/library")
# parameter_benchmark_archie.R — read-count permutation null model
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

source("scripts/parameter_benchmark_helpers.R")

NULL_NAME <- "archie_scramble"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dat <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl

scramble_counts <- function(dat, seed=42) {
  set.seed(seed)
  dat_scr <- dat
  idx <- sample(length(dat))
  mcols(dat_scr)$readsM <- mcols(dat)$readsM[idx]
  mcols(dat_scr)$readsN <- mcols(dat)$readsN[idx]
  message("  Original  mean methylation: ",
          round(sum(mcols(dat)$readsM)/sum(mcols(dat)$readsN)*100, 1), "%")
  message("  Scrambled mean methylation: ",
          round(sum(mcols(dat_scr)$readsM)/sum(mcols(dat_scr)$readsN)*100, 1), "%")
  dat_scr
}

message("Creating read-count permutation scrambled datasets...")
message("ASO_VPA:")
aso_vpa_scr  <- scramble_counts(aso_vpa,  seed=42)
message("ASO_CTRL:")
aso_ctrl_scr <- scramble_counts(aso_ctrl, seed=123)
message("  Done")

run_benchmark_loop(aso_vpa, aso_ctrl, aso_vpa_scr, aso_ctrl_scr, NULL_NAME)
