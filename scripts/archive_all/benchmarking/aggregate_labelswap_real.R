.libPaths("~/R/library")
source("scripts/benchmarking/labelswap_real_helpers.R")

# ---------------------------------------------------------------------------
# Aggregate per-task RDS files from the corrected label-swap array into:
#
# 1) parameter_benchmark_labelswap_real_<method>_ws<WS>_strict.csv
#    One row per (method, ws). Column format matches the ORIGINAL broken
#    parameter_benchmark_labelswap_*_ws*_strict.csv files exactly, so the
#    existing plot_radu_panels.R script will read these without modification.
#    n_scrambled = MEAN of perms 1..9. ratio = n_real / mean_scrambled.
#
# 2) labelswap_real_per_perm.csv
#    One row per (method, ws, perm_id) -- the full underlying distribution,
#    so you can show real value + scrambled spread (boxplot/violin/SD) on
#    the new panels.
#
# 3) labelswap_real_summary.csv
#    Compact summary (method, ws, n_real, n_scr_mean, n_scr_sd, n_scr_min,
#    n_scr_max, ratio, delta).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

files <- list.files(PERM_DIR, pattern = "^labelswap_real_.*\\.rds$",
                    full.names = TRUE)
message("Found ", length(files), " per-task RDS files in ", PERM_DIR)

if (length(files) == 0) stop("No RDS files to aggregate. Did the array run?")

raw <- do.call(rbind, lapply(files, readRDS))

# ---- per-perm CSV (one row per task) --------------------------------------
write.csv(raw, file.path(OUT_DIR, "labelswap_real_per_perm.csv"),
          row.names = FALSE)
message("Wrote: ", file.path(OUT_DIR, "labelswap_real_per_perm.csv"),
        " (", nrow(raw), " rows)")

# ---- check completeness ---------------------------------------------------
methods      <- c("bins", "neighbourhood", "noise_filter")
window_sizes <- c(100, 200, 300, 500, 1000, 2000)
expected <- expand.grid(method = methods, window_size = window_sizes,
                        perm_id = 0:9, stringsAsFactors = FALSE)
expected <- expected[!(expected$method == "neighbourhood" &
                       expected$window_size != window_sizes[1]), ]
missing <- anti_join(expected, raw, by = c("method","window_size","perm_id"))
if (nrow(missing) > 0) {
  message("WARNING: ", nrow(missing), " task(s) missing:")
  print(missing)
  writeLines(capture.output(print(missing, row.names = FALSE)),
             file.path(OUT_DIR, "missing_tasks.txt"))
} else {
  message("All ", nrow(expected), " expected tasks present.")
}

# ---- compact summary ------------------------------------------------------
real_rows <- raw %>% filter(is_real)
scr_rows  <- raw %>% filter(!is_real)

real_summary <- real_rows %>%
  select(method, window_size,
         n_real      = n_dmrs,
         size_real_median = size_median,
         size_real_mean   = size_mean,
         size_real_q1     = size_q1,
         size_real_q3     = size_q3,
         size_real_min    = size_min,
         size_real_max    = size_max)

scr_summary <- scr_rows %>%
  group_by(method, window_size) %>%
  summarise(
    n_perms          = n(),
    n_scr_mean       = mean(n_dmrs, na.rm = TRUE),
    n_scr_sd         = sd(n_dmrs, na.rm = TRUE),
    n_scr_min        = min(n_dmrs, na.rm = TRUE),
    n_scr_max        = max(n_dmrs, na.rm = TRUE),
    size_scr_median  = mean(size_median, na.rm = TRUE),
    size_scr_mean    = mean(size_mean,   na.rm = TRUE),
    size_scr_q1      = mean(size_q1,     na.rm = TRUE),
    size_scr_q3      = mean(size_q3,     na.rm = TRUE),
    size_scr_min     = mean(size_min,    na.rm = TRUE),
    size_scr_max     = mean(size_max,    na.rm = TRUE),
    .groups = "drop"
  )

compact <- real_summary %>%
  left_join(scr_summary, by = c("method","window_size")) %>%
  mutate(
    ratio = ifelse(n_scr_mean > 0, round(n_real / n_scr_mean, 3),
            ifelse(n_real > 0, Inf, NA_real_)),
    delta = n_real - n_scr_mean
  ) %>%
  arrange(method, window_size)

write.csv(compact, file.path(OUT_DIR, "labelswap_real_summary.csv"),
          row.names = FALSE)
message("Wrote: ", file.path(OUT_DIR, "labelswap_real_summary.csv"))
print(compact %>% select(method, window_size, n_real, n_scr_mean,
                         n_scr_sd, ratio, delta),
      n = 50)

# ---- per-cell CSVs matching original parameter_benchmark_labelswap format -
# Column schema (from make_result_row in parameter_benchmark_helpers.R):
#   method, window_size, mode, kernel, scramble_method,
#   n_real, n_scrambled, ratio,
#   size_real_*, size_scr_*,
#   t_real_secs, t_scr_secs, t_wall_secs,
#   error_real, error_scr
#
# We use MEAN of perms 1..9 as n_scrambled so the original plot_radu_panels.R
# reads it correctly. The full distribution lives in labelswap_real_per_perm.csv
# for any plot that wants to show spread.

OUT_LEGACY <- "results/dmr_benchmark"
dir.create(OUT_LEGACY, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(compact))) {
  row <- compact[i, ]
  m <- row$method;  w <- row$window_size
  ker_lab <- if (m == "noise_filter") "triangular" else "NA"

  # Per-cell timings: average across perm_ids for that cell
  cell_raw <- raw %>% filter(method == m, window_size == w)
  t_real <- cell_raw %>% filter(is_real)  %>% pull(t_dmr_secs)
  t_scr  <- cell_raw %>% filter(!is_real) %>% pull(t_dmr_secs)

  legacy <- data.frame(
    method            = m,
    window_size       = w,
    mode              = "strict",
    kernel            = ker_lab,
    scramble_method   = "label_swap_real",
    n_real            = row$n_real,
    n_scrambled       = round(row$n_scr_mean, 1),
    ratio             = row$ratio,
    size_real_median  = row$size_real_median,
    size_real_mean    = row$size_real_mean,
    size_real_q1      = row$size_real_q1,
    size_real_q3      = row$size_real_q3,
    size_real_min     = row$size_real_min,
    size_real_max     = row$size_real_max,
    size_scr_median   = round(row$size_scr_median, 1),
    size_scr_mean     = round(row$size_scr_mean,   1),
    size_scr_q1       = round(row$size_scr_q1,     1),
    size_scr_q3       = round(row$size_scr_q3,     1),
    size_scr_min      = round(row$size_scr_min,    1),
    size_scr_max      = round(row$size_scr_max,    1),
    t_real_secs       = if (length(t_real) > 0) round(t_real[1], 1) else NA_real_,
    t_scr_secs        = if (length(t_scr)  > 0) round(mean(t_scr, na.rm = TRUE), 1) else NA_real_,
    t_wall_secs       = NA_real_,
    error_real        = NA_character_,
    error_scr         = NA_character_,
    n_scr_sd          = round(row$n_scr_sd, 1),
    n_perms           = row$n_perms,
    stringsAsFactors  = FALSE
  )

  out_path <- file.path(OUT_LEGACY,
    sprintf("parameter_benchmark_labelswap_real_%s_ws%d_strict.csv", m, w))
  write.csv(legacy, out_path, row.names = FALSE)
  message("Wrote legacy-format CSV: ", basename(out_path))
}

message("\nDone. Three output sets:")
message("  ", OUT_DIR,   "/labelswap_real_per_perm.csv  (one row per task)")
message("  ", OUT_DIR,   "/labelswap_real_summary.csv   (compact summary)")
message("  ", OUT_LEGACY,"/parameter_benchmark_labelswap_real_*_strict.csv  (legacy-schema, one per cell)")
