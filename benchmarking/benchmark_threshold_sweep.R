.libPaths("~/R/library")
source("scripts/benchmarking/parameter_benchmark_helpers.R")

# Threshold sweep: minProportionDifference = 0.1, 0.2, 0.3, 0.4
# Requested by Radu for the label-swap and read-count permutation panels,
# y-axis = number of DMRs (not coverage). Reuses run_dmrs_one(),
# run_dmrs_pair_parallel(), make_result_row(), load_chr1_data() exactly as
# the original strict/loose scripts do -- only minDiff is overridden.
#
# Scope: bins method only (matches what plot_radu_panels.R actually plots;
# neighbourhood and noise_filter already shown separately/excluded per
# earlier findings), 6 window sizes, 4 thresholds, 2 null models.
# 6 x 4 x 2 = 48 tasks.

run_dmrs_one_threshold <- function(treat, ctrl, ws, min_diff, region) {
  th <- get_thresholds("bins", strict = TRUE, ws = ws)  # base: strict pval/minCpG/minGap
  th$minDiff <- min_diff                                  # override only this

  setTimeLimit(cpu = 1800, elapsed = 1800, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE))
  tryCatch({
    computeDMRs(treat, ctrl, regions = region, context = "CG",
      method = "bins", binSize = ws, test = th$test, pValueThreshold = th$pval,
      minCytosinesCount = th$minCpG, minProportionDifference = th$minDiff,
      minGap = th$minGap, minSize = th$minSize, minReadsPerCytosine = 4,
      parallel = FALSE)
  }, error = function(e) {
    message("  Error: ", e$message)
    structure(list(), class = "dmr_error", message = e$message)
  })
}

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

window_sizes <- c(100, 200, 300, 500, 1000, 2000)
thresholds   <- c(0.1, 0.2, 0.3, 0.4)
null_models  <- c("label_swap", "read_count_permutation")

combos <- expand.grid(window_size = window_sizes, min_diff = thresholds,
                      null_model = null_models, stringsAsFactors = FALSE)
message("Total tasks: ", nrow(combos))

if (TASK_ID > nrow(combos)) {
  message("TASK_ID ", TASK_ID, " exceeds combo table (nrow=", nrow(combos), ") -- skipping")
  quit(status = 0)
}

ws         <- combos$window_size[TASK_ID]
min_diff   <- combos$min_diff[TASK_ID]
null_model <- combos$null_model[TASK_ID]
message("Task ", TASK_ID, ": ws=", ws, " minDiff=", min_diff, " null=", null_model)

dat      <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl
region   <- GRanges(CHROM, IRanges(1, REGION_END))

OUT_DIR_SWEEP <- "results/dmr_benchmark_threshold_sweep"
dir.create(OUT_DIR_SWEEP, recursive = TRUE, showWarnings = FALSE)

if (null_model == "label_swap") {
  aso_vpa_scr  <- aso_ctrl
  aso_ctrl_scr <- aso_vpa

  t_start <- Sys.time()
  out_real <- run_dmrs_one_threshold(aso_vpa, aso_ctrl, ws, min_diff, region)
  out_scr  <- run_dmrs_one_threshold(aso_vpa_scr, aso_ctrl_scr, ws, min_diff, region)
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

} else {
  scramble_data <- function(dat1, dat2, seed = 42) {
    set.seed(seed)
    idx1 <- sample(seq_len(length(dat1)), replace = FALSE)
    idx2 <- sample(seq_len(length(dat2)), replace = FALSE)
    dat1_scr <- dat1; dat2_scr <- dat2
    mcols(dat1_scr)$readsM <- mcols(dat1)$readsM[idx1]
    mcols(dat1_scr)$readsN <- mcols(dat1)$readsN[idx1]
    mcols(dat2_scr)$readsM <- mcols(dat2)$readsM[idx2]
    mcols(dat2_scr)$readsN <- mcols(dat2)$readsN[idx2]
    list(dat1 = dat1_scr, dat2 = dat2_scr)
  }

  N_PERMS <- 20
  t_start <- Sys.time()
  out_real <- run_dmrs_one_threshold(aso_vpa, aso_ctrl, ws, min_diff, region)
  n_real <- if (is_error(out_real)) NA_integer_ else length(out_real)
  message("  Real DMRs: ", n_real)

  scr_counts <- integer(N_PERMS)
  for (s in 1:N_PERMS) {
    scr <- scramble_data(aso_vpa, aso_ctrl, seed = s)
    d_s <- run_dmrs_one_threshold(scr$dat1, scr$dat2, ws, min_diff, region)
    scr_counts[s] <- if (is_error(d_s)) NA_integer_ else length(d_s)
    message("  seed ", s, "/", N_PERMS, ": ", scr_counts[s])
  }
  t_wall <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  mean_scr <- mean(scr_counts, na.rm = TRUE)
  sd_scr   <- sd(scr_counts, na.rm = TRUE)
  ratio    <- if (mean_scr == 0) Inf else round(n_real / mean_scr, 3)

  row <- data.frame(
    method = "bins", window_size = ws, mode = paste0("minDiff_", min_diff),
    minDiff = min_diff, n_real = n_real,
    mean_scrambled = mean_scr, sd_scrambled = sd_scr, ratio = ratio,
    n_seeds = N_PERMS, scramble_method = "read_count_permutation",
    t_wall_secs = round(t_wall, 1), stringsAsFactors = FALSE
  )
  out_csv <- file.path(OUT_DIR_SWEEP,
    sprintf("sweep_readcount_bins_ws%d_md%s.csv", ws, gsub("\\.", "", as.character(min_diff))))
  write.csv(row, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
}
