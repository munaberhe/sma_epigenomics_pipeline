# parameter_benchmark_turbo.R
# PLAN C: label-swap scramble, 32-worker serial real+scrambled.
source("scripts/parameter_benchmark_helpers_turbo.R")
set.seed(42)

dat <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl

COV_DIR <- "results/alignments/bs/by_chr"
CHROM   <- "chr1"

load_one <- function(samples) {
  paths <- file.path(COV_DIR, paste0(samples, "_", CHROM, ".CpG_report.txt.gz"))
  paths <- paths[file.exists(paths)]
  readBismarkPool(paths)
}

aso_vpa_scr  <- load_one(c("ASO_CTRL_1", "ASO_VPA_2",  "ASO_VPA_3"))
aso_ctrl_scr <- load_one(c("ASO_VPA_1",  "ASO_CTRL_2", "ASO_CTRL_3"))

run_benchmark_loop(aso_vpa, aso_ctrl, aso_vpa_scr, aso_ctrl_scr, null_name = "label_swap")
