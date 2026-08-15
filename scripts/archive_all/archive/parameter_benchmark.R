.libPaths("~/R/library")
# parameter_benchmark.R — label-swap null model
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

source("scripts/parameter_benchmark_helpers.R")

NULL_NAME <- "label_swap"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dat <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl

message("Creating label-swap scrambled datasets...")
aso_vpa_scr  <- aso_ctrl
aso_ctrl_scr <- aso_vpa
message("  Done")

run_benchmark_loop(aso_vpa, aso_ctrl, aso_vpa_scr, aso_ctrl_scr, NULL_NAME)
