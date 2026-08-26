.libPaths("~/R/library")
source("scripts/benchmarking/parameter_benchmark_helpers.R")

# Targeted rerun of ONLY the combinations that failed in
# benchmark_threshold_sweep.R due to the 1800s internal time limit being
# too short at parallel=FALSE. This version enables parallel=TRUE, cores=32
# (matching the original benchmark_readcount_array.R convention) and raises
# the internal time limit, since CPU parallelism should substantially speed
# up computeDMRs() for larger window sizes.
#
# Failed combinations from the first sweep (label_swap only; read_count
# permutation is not attempted here -- see separate note):
#   ws=100,  md=0.2
#   ws=200,  md=0.2
#   ws=300,  md=0.2
#   ws=500,  md=0.2, 0.3, 0.4
#   ws=1000, md=0.2, 0.3, 0.4
#   ws=2000, md=0.2, 0.3, 0.4
# 13 failed combinations total -- this array covers exactly those 13 tasks.

FAILED_COMBOS <- data.frame(
  window_size = c(100, 200, 300, 500, 500, 500, 1000, 1000, 1000, 2000, 2000, 2000, 500),
  min_diff    = c(0.2, 0.2, 0.2, 0.2, 0.3, 0.4, 0.2,  0.3,  0.4,  0.2,  0.3,  0.4,  0.3)
)
FAILED_COMBOS <- unique(FAILED_COMBOS)
message("Total tasks: ", nrow(FAILED_COMBOS))

run_dmrs_one_parallel <- function(treat, ctrl, ws, min_diff, region, cores = 16) {
  th <- get_thresholds("bins", strict = TRUE, ws = ws)
  th$minDiff <- min_diff

  # Raised from 1800s -- parallel processing should make this feasible,
  # but give real headroom since larger window sizes are still slower.
  setTimeLimit(cpu = 18000, elapsed = 18000, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE))
  tryCatch({
    computeDMRs(treat, ctrl, regions = region, context = "CG",
      method = "bins", binSize = ws, test = th$test, pValueThreshold = th$pval,
      minCytosinesCount = th$minCpG, minProportionDifference = th$minDiff,
      minGap = th$minGap, minSize = th$minSize, minReadsPerCytosine = 4,
      parallel = TRUE, cores = cores)
  }, error = function(e) {
    message("  Error: ", e$message)
    structure(list(), class = "dmr_error", message = e$message)
  })
}

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
if (TASK_ID > nrow(FAILED_COMBOS)) {
  message("TASK_ID ", TASK_ID, " exceeds combo table -- skipping")
  quit(status = 0)
}

ws       <- FAILED_COMBOS$window_size[TASK_ID]
min_diff <- FAILED_COMBOS$min_diff[TASK_ID]
message("Task ", TASK_ID, ": ws=", ws, " minDiff=", min_diff, " (parallel, cores=32)")

dat      <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl
region   <- GRanges(CHROM, IRanges(1, REGION_END))

aso_vpa_scr  <- aso_ctrl
aso_ctrl_scr <- aso_vpa

OUT_DIR_SWEEP <- "results/dmr_benchmark_threshold_sweep"
dir.create(OUT_DIR_SWEEP, recursive = TRUE, showWarnings = FALSE)

t_start <- Sys.time()
out_real <- run_dmrs_one_parallel(aso_vpa, aso_ctrl, ws, min_diff, region)
out_scr  <- run_dmrs_one_parallel(aso_vpa_scr, aso_ctrl_scr, ws, min_diff, region)
t_wall <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

n_real <- if (is_error(out_real)) NA_integer_ else length(out_real)
n_scr  <- if (is_error(out_scr))  NA_integer_ else length(out_scr)
message("  Real: ", n_real, "  Scrambled: ", n_scr, "  wall=", round(t_wall), "s")

pair <- list(real = out_real, scr = out_scr, t_real = t_wall, t_scr = t_wall,
            t_wall = t_wall)
row <- make_result_row(pair, "bins", ws, paste0("minDiff_", min_diff), "NA", "label_swap")
row$minDiff <- min_diff
out_csv <- file.path(OUT_DIR_SWEEP,
  sprintf("sweep_labelswap_bins_ws%d_md%s.csv", ws, gsub("\\.", "", as.character(min_diff))))
write.csv(row, out_csv, row.names = FALSE)
message("Saved: ", out_csv)
