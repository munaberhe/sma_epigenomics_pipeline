.libPaths("~/R/library")
# parameter_benchmark_serial.R — PLAN B (label-swap)
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

.this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) > 0) dirname(normalizePath(f)) else getwd()
  }
)
source(file.path(.this_dir, "parameter_benchmark_helpers_serial.R"))

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
