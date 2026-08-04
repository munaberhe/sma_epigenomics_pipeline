.libPaths("~/R/library")
source("scripts/benchmarking/parameter_benchmark_helpers.R")

# ---------------------------------------------------------------------------
# Intermediate-threshold (minProportionDifference = 0.02) benchmark, matching
# Radu's selection from the SMN2 sensitivity scan. Reuses every existing,
# validated function from parameter_benchmark_helpers.R unchanged --
# run_dmrs_one(), run_dmrs_pair_parallel(), make_result_row() are called
# exactly as the original label-swap and read-count scripts call them.
#
# The ONLY new piece is a thin wrapper around run_dmrs_one() that overrides
# minProportionDifference to 0.02 after get_thresholds() returns its
# strict-mode defaults (pval=0.01, minCpG=4, minGap=200, minSize as before).
# This keeps get_thresholds() itself, and every other script that calls it,
# completely unchanged -- existing strict/loose results are not at risk.
#
# Restricted to bins and neighbourhood (noise_filter excluded: confirmed
# degenerate -- flat at zero -- under both strict and loose thresholds in
# the existing benchmark CSVs, so not expected to behave differently here).
#
# Both null models included as separate array ranges so this one script
# covers what was previously two separate scripts:
#   tasks 1-12:  label_swap      (6 window sizes x 2 methods)
#   tasks 13-24: read_count_perm (6 window sizes x 2 methods, 20 perms each,
#                matching benchmark_readcount_array.R's averaging convention)
# ---------------------------------------------------------------------------

run_dmrs_one_intermediate <- function(treat, ctrl, method, ws, region) {
  th <- get_thresholds(method, strict = TRUE, ws = ws)  # base: strict mode
  th$minDiff <- 0.02                                      # override only this

  setTimeLimit(cpu = 1800, elapsed = 1800, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE))
  tryCatch({
    if (method == "bins") {
      computeDMRs(treat, ctrl, regions = region, context = "CG",
        method = "bins", binSize = ws, test = th$test, pValueThreshold = th$pval,
        minCytosinesCount = th$minCpG, minProportionDifference = th$minDiff,
        minGap = th$minGap, minSize = th$minSize, minReadsPerCytosine = 4,
        parallel = FALSE)
    } else {
      computeDMRs(treat, ctrl, regions = region, context = "CG",
        method = "neighbourhood", test = th$test, pValueThreshold = th$pval,
        minCytosinesCount = th$minCpG, minProportionDifference = th$minDiff,
        minGap = th$minGap, minSize = th$minSize, minReadsPerCytosine = 4,
        parallel = FALSE)
    }
  }, error = function(e) {
    message("  Error: ", e$message)
    structure(list(), class = "dmr_error", message = e$message)
  })
}

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

methods      <- c("bins", "neighbourhood")
window_sizes <- c(100, 200, 300, 500, 1000, 2000)

ls_combos <- expand.grid(method = methods, window_size = window_sizes,
                         null_model = "label_swap", stringsAsFactors = FALSE)
ls_combos <- ls_combos[!(ls_combos$method == "neighbourhood" &
                         ls_combos$window_size != window_sizes[1]), ]

rc_combos <- expand.grid(method = methods, window_size = window_sizes,
                         null_model = "read_count_permutation", stringsAsFactors = FALSE)
rc_combos <- rc_combos[!(rc_combos$method == "neighbourhood" &
                         rc_combos$window_size != window_sizes[1]), ]

combos <- rbind(ls_combos, rc_combos)
message("Total tasks: ", nrow(combos))

if (TASK_ID > nrow(combos)) {
  message("TASK_ID ", TASK_ID, " exceeds combo table (nrow=", nrow(combos), ") -- skipping")
  quit(status = 0)
}

method     <- combos$method[TASK_ID]
ws         <- combos$window_size[TASK_ID]
null_model <- combos$null_model[TASK_ID]
message("Task ", TASK_ID, ": ", method, " ws=", ws, " null=", null_model,
        " minDiff=0.02")

dat      <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl
region   <- GRanges(CHROM, IRanges(1, REGION_END))

OUT_DIR_2PCT <- "results/dmr_benchmark_2pct"
dir.create(OUT_DIR_2PCT, recursive = TRUE, showWarnings = FALSE)

if (null_model == "label_swap") {
  # exact same construction as benchmark_labelswap_array.R: swap which
  # condition is treatment vs control
  aso_vpa_scr  <- aso_ctrl
  aso_ctrl_scr <- aso_vpa

  t_start <- Sys.time()
  out_real <- run_dmrs_one_intermediate(aso_vpa, aso_ctrl, method, ws, region)
  out_scr  <- run_dmrs_one_intermediate(aso_vpa_scr, aso_ctrl_scr, method, ws, region)
  t_wall <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  n_real <- if (is_error(out_real)) NA_integer_ else length(out_real)
  n_scr  <- if (is_error(out_scr))  NA_integer_ else length(out_scr)
  message("  Real: ", n_real, "  Scrambled: ", n_scr, "  wall=", round(t_wall), "s")

  pair <- list(real = out_real, scr = out_scr, t_real = t_wall, t_scr = t_wall,
              t_wall = t_wall)
  row <- make_result_row(pair, method, ws, "intermediate_2pct", "NA", "label_swap")
  out_csv <- file.path(OUT_DIR_2PCT,
    sprintf("benchmark_2pct_labelswap_%s_ws%d.csv", method, ws))
  write.csv(row, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)

} else {
  # exact same construction as benchmark_readcount_array.R: shuffle
  # readsM/readsN pairing within each condition independently, 20 seeds,
  # averaged
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
  out_real <- run_dmrs_one_intermediate(aso_vpa, aso_ctrl, method, ws, region)
  n_real <- if (is_error(out_real)) NA_integer_ else length(out_real)
  message("  Real DMRs: ", n_real)

  scr_counts <- integer(N_PERMS)
  for (s in 1:N_PERMS) {
    scr <- scramble_data(aso_vpa, aso_ctrl, seed = s)
    d_s <- run_dmrs_one_intermediate(scr$dat1, scr$dat2, method, ws, region)
    scr_counts[s] <- if (is_error(d_s)) NA_integer_ else length(d_s)
    message("  seed ", s, "/", N_PERMS, ": ", scr_counts[s])
  }
  t_wall <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  mean_scr <- mean(scr_counts, na.rm = TRUE)
  sd_scr   <- sd(scr_counts, na.rm = TRUE)
  ratio    <- if (mean_scr == 0) Inf else round(n_real / mean_scr, 3)

  row <- data.frame(
    method = method, window_size = ws, mode = "intermediate_2pct",
    minDiff = 0.02, n_real = n_real,
    mean_scrambled = mean_scr, sd_scrambled = sd_scr, ratio = ratio,
    ratio_lower = round(n_real / (mean_scr + sd_scr), 3),
    ratio_upper = round(n_real / pmax(mean_scr - sd_scr, 0.001), 3),
    n_seeds = N_PERMS, scramble_method = "read_count_permutation",
    t_wall_secs = round(t_wall, 1), stringsAsFactors = FALSE
  )
  out_csv <- file.path(OUT_DIR_2PCT,
    sprintf("benchmark_2pct_readcount_%s_ws%d.csv", method, ws))
  write.csv(row, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
}
