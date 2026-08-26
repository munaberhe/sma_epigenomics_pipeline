.libPaths("~/R/library")
source("scripts/benchmarking/labelswap_real_helpers.R")

# Per-task script for the corrected label-swap array.
# Each task = one (method, window, perm_id) cell.
# Total: 3 methods x 6 windows x 10 perms = 180 tasks.
# Note: neighbourhood is window-invariant (DMRcaller's neighbourhood method
# does not take a windowSize argument), so for neighbourhood we only need
# ws = 100 (the first window). We mark the remaining neighbourhood tasks
# as skipped via early exit, matching the original benchmark_labelswap_array.R
# logic at line 17.
#
# After all tasks finish, run aggregate_labelswap_real.R to assemble the
# CSV summary in the same column format as the original
# parameter_benchmark_labelswap_*_strict.csv files.

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

methods      <- c("bins", "neighbourhood", "noise_filter")
window_sizes <- c(100, 200, 300, 500, 1000, 2000)
n_perms      <- 10  # perm_id 0..9

combos <- expand.grid(
  method      = methods,
  window_size = window_sizes,
  perm_id     = 0:(n_perms - 1),
  stringsAsFactors = FALSE
)

# Drop neighbourhood tasks for ws != 100 (window-invariant, matches original)
combos <- combos[!(combos$method == "neighbourhood" &
                   combos$window_size != window_sizes[1]), ]
# Total now: 6*10 + 1*10 + 6*10 = 130 tasks. Submit array as 1-130.

if (TASK_ID > nrow(combos)) {
  message("TASK_ID ", TASK_ID, " exceeds combo table (nrow=", nrow(combos),
          ") -- skipping")
  quit(status = 0)
}

method  <- combos$method[TASK_ID]
ws      <- combos$window_size[TASK_ID]
perm_id <- combos$perm_id[TASK_ID]

# Output path: one RDS per task. Aggregator scans these.
out_rds <- file.path(PERM_DIR,
  sprintf("labelswap_real_%s_ws%d_perm%02d.rds", method, ws, perm_id))

if (file.exists(out_rds)) {
  message("Already done: ", basename(out_rds), " -- skipping")
  quit(status = 0)
}

message("Task ", TASK_ID, ": method=", method, " ws=", ws,
        " perm_id=", perm_id)

# Enumerate partitions and pick this perm's group assignment
partitions <- enumerate_partitions()
perm_row   <- partitions[partitions$perm_id == perm_id, ]
stopifnot(nrow(perm_row) == 1)
message("  groupA: ", perm_row$groupA)
message("  groupB: ", perm_row$groupB)

# Pool the replicates for this partition
t0 <- Sys.time()
pair <- build_pair_for_perm(perm_row)
t_pool <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message("  Pool wallclock: ", round(t_pool), "s   |  treat CpGs=",
        length(pair$treat), "  ctrl CpGs=", length(pair$ctrl))

# Run one DMR call
region <- GRanges(CHROM, IRanges(1, REGION_END))
ker    <- if (method == "noise_filter") "triangular" else "NA"

t1 <- Sys.time()
dmrs <- run_dmr_call(pair$treat, pair$ctrl, method, ws, ker, region)
t_dmr <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

n_dmrs <- if (is_error(dmrs)) NA_integer_ else length(dmrs)
err    <- if (is_error(dmrs)) attr(dmrs, "message") else NA_character_
sz     <- dmr_size_summary(dmrs)
message("  DMRs: ", n_dmrs, "   DMR call wallclock: ", round(t_dmr), "s")

res <- data.frame(
  method      = method,
  window_size = ws,
  perm_id     = perm_id,
  is_real     = perm_row$is_real,
  groupA      = perm_row$groupA,
  groupB      = perm_row$groupB,
  kernel      = ker,
  n_dmrs      = n_dmrs,
  size_median = sz$median,
  size_mean   = sz$mean,
  size_q1     = sz$q1,
  size_q3     = sz$q3,
  size_min    = sz$min,
  size_max    = sz$max,
  t_pool_secs = round(t_pool, 1),
  t_dmr_secs  = round(t_dmr, 1),
  error       = err,
  stringsAsFactors = FALSE
)

saveRDS(res, out_rds)
message("Saved: ", out_rds)
