.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(dplyr)
})

# Aggregate per-task RDS files from the fine-grained read-count permutation
# array into a compact summary, matching the structure of
# labelswap_real_summary.csv / aggregate_labelswap_real.R.
#
# Each task RDS has: window_size, minDiff, perm_id, is_real, n_dmrs,
# size_*, t_wall_secs, error.
# perm_id == 0 is the real (unpermuted) call; perm_id 1-20 are the 20
# read-count-permuted nulls.

PERM_DIR <- "results/dmr_benchmark_readcount_real/perm_rds"
OUT_DIR  <- "results/dmr_benchmark_readcount_real"

files <- list.files(PERM_DIR, pattern = "^readcount_real_.*\\.rds$", full.names = TRUE)
message("Found ", length(files), " per-task RDS files in ", PERM_DIR)

if (length(files) == 0) stop("No RDS files to aggregate. Did the array run?")

raw <- do.call(rbind, lapply(files, readRDS))
write.csv(raw, file.path(OUT_DIR, "readcount_real_per_perm.csv"), row.names = FALSE)
message("Wrote: ", file.path(OUT_DIR, "readcount_real_per_perm.csv"), " (", nrow(raw), " rows)")

# check completeness
window_sizes <- c(100, 200, 300, 500, 1000, 2000)
thresholds   <- c(0.1, 0.2, 0.3, 0.4)
perm_ids     <- 0:20
expected <- expand.grid(window_size = window_sizes, minDiff = thresholds,
                        perm_id = perm_ids, stringsAsFactors = FALSE)
missing <- anti_join(expected, raw, by = c("window_size", "minDiff", "perm_id"))
if (nrow(missing) > 0) {
  message("WARNING: ", nrow(missing), " task(s) missing:")
  print(missing)
  writeLines(capture.output(print(missing, row.names = FALSE)),
             file.path(OUT_DIR, "missing_tasks.txt"))
} else {
  message("All ", nrow(expected), " expected tasks present.")
}

# compact summary
real_rows <- raw %>% filter(is_real)
scr_rows  <- raw %>% filter(!is_real)

real_summary <- real_rows %>%
  select(window_size, minDiff,
         n_real = n_dmrs,
         size_real_median = size_median, size_real_mean = size_mean,
         size_real_q1 = size_q1, size_real_q3 = size_q3,
         size_real_min = size_min, size_real_max = size_max)

scr_summary <- scr_rows %>%
  group_by(window_size, minDiff) %>%
  summarise(
    n_perms    = n(),
    n_scr_mean = mean(n_dmrs, na.rm = TRUE),
    n_scr_sd   = sd(n_dmrs, na.rm = TRUE),
    n_scr_min  = suppressWarnings(min(n_dmrs, na.rm = TRUE)),
    n_scr_max  = suppressWarnings(max(n_dmrs, na.rm = TRUE)),
    .groups = "drop"
  )

compact <- real_summary %>%
  left_join(scr_summary, by = c("window_size", "minDiff")) %>%
  mutate(
    ratio = ifelse(!is.na(n_scr_mean) & n_scr_mean > 0, round(n_real / n_scr_mean, 3),
            ifelse(!is.na(n_real) & n_real > 0, Inf, NA_real_)),
    delta = n_real - n_scr_mean
  ) %>%
  arrange(window_size, minDiff)

write.csv(compact, file.path(OUT_DIR, "readcount_real_summary.csv"), row.names = FALSE)
message("Wrote: ", file.path(OUT_DIR, "readcount_real_summary.csv"))

cat("\n=== Summary ===\n")
print(as.data.frame(compact[, c("window_size","minDiff","n_real","n_scr_mean","n_scr_sd","ratio","delta")]),
      row.names = FALSE)

# legacy-format CSVs (one per window/threshold cell), matching the
# schema plot_radu_panels.R / the benchmark CSVs use
OUT_LEGACY <- "results/dmr_benchmark_threshold_sweep"
dir.create(OUT_LEGACY, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(compact))) {
  row <- compact[i, ]
  w <- row$window_size; md <- row$minDiff
  legacy <- data.frame(
    method = "bins", window_size = w, mode = paste0("minDiff_", md), minDiff = md,
    n_real = row$n_real, mean_scrambled = round(row$n_scr_mean, 1),
    sd_scrambled = round(row$n_scr_sd, 1), ratio = row$ratio,
    n_seeds = row$n_perms, scramble_method = "read_count_permutation",
    stringsAsFactors = FALSE
  )
  out_path <- file.path(OUT_LEGACY,
    sprintf("sweep_readcount_real_bins_ws%d_md%s.csv", w, gsub("\\.", "", as.character(md))))
  write.csv(legacy, out_path, row.names = FALSE)
  message("Wrote: ", basename(out_path))
}

message("\nDone.")
