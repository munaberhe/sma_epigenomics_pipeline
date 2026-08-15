.libPaths("~/R/library")
source("scripts/benchmarking/parameter_benchmark_helpers.R")

PERM_DIR <- "results/dmr_benchmark_readcount_real/perm_rds"
dir.create(PERM_DIR, recursive = TRUE, showWarnings = FALSE)

run_dmrs_one_parallel <- function(treat, ctrl, ws, min_diff, region, cores = 16) {
  th <- get_thresholds("bins", strict = TRUE, ws = ws)
  th$minDiff <- min_diff
  setTimeLimit(cpu = 7200, elapsed = 7200, transient = TRUE)
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

scramble_data <- function(dat1, dat2, seed) {
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

dmr_size_summary <- function(dmrs) {
  if (is_error(dmrs) || length(dmrs) == 0)
    return(list(median = NA_real_, mean = NA_real_, q1 = NA_real_,
                q3 = NA_real_, min = NA_real_, max = NA_real_))
  w <- as.numeric(width(dmrs))
  list(median = median(w), mean = round(mean(w), 1),
       q1 = as.numeric(quantile(w, 0.25)), q3 = as.numeric(quantile(w, 0.75)),
       min = min(w), max = max(w))
}

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

window_sizes <- c(100, 200, 300, 500, 1000, 2000)
thresholds   <- c(0.1, 0.2, 0.3, 0.4)
perm_ids     <- 0:20

combos <- expand.grid(window_size = window_sizes, min_diff = thresholds,
                      perm_id = perm_ids, stringsAsFactors = FALSE)
message("Total tasks: ", nrow(combos))

if (TASK_ID > nrow(combos)) {
  message("TASK_ID ", TASK_ID, " exceeds combo table -- skipping")
  quit(status = 0)
}

ws       <- combos$window_size[TASK_ID]
min_diff <- combos$min_diff[TASK_ID]
perm_id  <- combos$perm_id[TASK_ID]
is_real  <- perm_id == 0
message("Task ", TASK_ID, ": ws=", ws, " minDiff=", min_diff,
        " perm_id=", perm_id, " (real=", is_real, ")")

out_rds <- file.path(PERM_DIR,
  sprintf("readcount_real_ws%d_md%s_perm%02d.rds",
          ws, gsub("\\.", "", as.character(min_diff)), perm_id))

if (file.exists(out_rds)) {
  message("Already done: ", basename(out_rds), " -- skipping")
  quit(status = 0)
}

dat      <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl
region   <- GRanges(CHROM, IRanges(1, REGION_END))

t0 <- Sys.time()
if (is_real) {
  dmrs <- run_dmrs_one_parallel(aso_vpa, aso_ctrl, ws, min_diff, region)
} else {
  scr <- scramble_data(aso_vpa, aso_ctrl, seed = perm_id)
  dmrs <- run_dmrs_one_parallel(scr$dat1, scr$dat2, ws, min_diff, region)
}
t_wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

n_dmrs <- if (is_error(dmrs)) NA_integer_ else length(dmrs)
err    <- if (is_error(dmrs)) attr(dmrs, "message") else NA_character_
sz     <- dmr_size_summary(dmrs)
message("  DMRs: ", n_dmrs, "  wall=", round(t_wall), "s")

res <- data.frame(
  window_size = ws, minDiff = min_diff, perm_id = perm_id, is_real = is_real,
  n_dmrs = n_dmrs, size_median = sz$median, size_mean = sz$mean,
  size_q1 = sz$q1, size_q3 = sz$q3, size_min = sz$min, size_max = sz$max,
  t_wall_secs = round(t_wall, 1), error = err, stringsAsFactors = FALSE
)
saveRDS(res, out_rds)
message("Saved: ", out_rds)
