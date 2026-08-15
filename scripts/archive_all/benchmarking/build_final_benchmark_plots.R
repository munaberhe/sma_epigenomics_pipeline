.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/dmr_benchmark_final_plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

NAVY <- "#1B4F8A"
RED  <- "#B2182B"
GREY <- "#6B7280"

# ---------------------------------------------------------------------------
# PLOT 1: corrected label-swap (real replicate-partition permutation)
# ---------------------------------------------------------------------------
ls_df <- read.csv("results/dmr_benchmark_labelswap_real/labelswap_real_summary.csv")
ls_df <- ls_df[ls_df$method == "bins", ]   # noise_filter excluded -- confirmed degenerate
ls_df <- ls_df[!is.na(ls_df$n_real), ]
ls_df$window_size <- as.numeric(ls_df$window_size)

ls_long <- rbind(
  data.frame(window_size = ls_df$window_size * 0.95, n = ls_df$n_real,
            series = "Real", stringsAsFactors = FALSE),
  data.frame(window_size = ls_df$window_size * 1.05, n = ls_df$n_scr_mean,
            series = "Scrambled (mean of 9 real partitions)", stringsAsFactors = FALSE)
)
ls_long$series <- factor(ls_long$series,
  levels = c("Real", "Scrambled (mean of 9 real partitions)"))

p_ls <- ggplot(ls_long, aes(x = window_size, y = n, colour = series, shape = series)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3, fill = "white", stroke = 0.9) +
  scale_colour_manual(values = c("Real" = NAVY, "Scrambled (mean of 9 real partitions)" = GREY)) +
  scale_shape_manual(values = c("Real" = 15, "Scrambled (mean of 9 real partitions)" = 0)) +
  scale_x_log10(breaks = c(100, 200, 300, 500, 1000, 2000)) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Label-swap: real replicate-partition permutation",
    subtitle = "DMRcaller-B, strict thresholds. Null = all 9 non-identity 3-vs-3 replicate\npartitions of the 6 ASO_VPA/ASO_CTRL replicates (exhaustive, not sampled).",
    x = "Bin/window size (bp)", y = "Number of DMRs (log scale)",
    colour = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "labelswap_real_corrected.pdf"), p_ls,
       width = 9, height = 6.5, device = cairo_pdf)
message("Saved: labelswap_real_corrected.pdf")

# ---------------------------------------------------------------------------
# PLOT 2: corrected read-count permutation, threshold sweep
# ---------------------------------------------------------------------------
rc_df <- read.csv("results/dmr_benchmark_readcount_real/readcount_real_summary.csv")
rc_df <- rc_df[!is.na(rc_df$n_real) & !is.na(rc_df$n_scr_mean), ]
rc_df$window_size <- as.numeric(rc_df$window_size)
rc_df$minDiff_lab <- paste0("minDiff = ", rc_df$minDiff)
rc_df$minDiff_lab <- factor(rc_df$minDiff_lab,
  levels = paste0("minDiff = ", c(0.1, 0.2, 0.3, 0.4)))

rc_long <- rbind(
  data.frame(window_size = rc_df$window_size * 0.95, n = rc_df$n_real,
            minDiff_lab = rc_df$minDiff_lab, series = "Real", stringsAsFactors = FALSE),
  data.frame(window_size = rc_df$window_size * 1.05, n = rc_df$n_scr_mean,
            minDiff_lab = rc_df$minDiff_lab,
            series = "Scrambled (mean of 20 read-count permutations)",
            stringsAsFactors = FALSE)
)
rc_long$series <- factor(rc_long$series,
  levels = c("Real", "Scrambled (mean of 20 read-count permutations)"))

p_rc <- ggplot(rc_long, aes(x = window_size, y = n, colour = series, shape = series)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.8, fill = "white", stroke = 0.8) +
  scale_colour_manual(values = c("Real" = RED, "Scrambled (mean of 20 read-count permutations)" = GREY)) +
  scale_shape_manual(values = c("Real" = 17, "Scrambled (mean of 20 read-count permutations)" = 2)) +
  scale_x_log10(breaks = c(100, 200, 300, 500, 1000, 2000)) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~minDiff_lab, nrow = 1) +
  labs(
    title = "Read-count permutation: threshold sweep",
    subtitle = "DMRcaller-B. Null = 20 permutations shuffling readsM/readsN pairing\nwithin each condition, independently, preserving per-CpG coverage.",
    x = "Bin/window size (bp)", y = "Number of DMRs (log scale)",
    colour = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(OUT_DIR, "readcount_real_corrected.pdf"), p_rc,
       width = 14, height = 5.5, device = cairo_pdf)
message("Saved: readcount_real_corrected.pdf")

message("\nBoth plots saved in: ", OUT_DIR)
