.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/dmr_benchmark_final_plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

NAVY <- "#1B4F8A"
RED  <- "#B2182B"
GREY <- "#6B7280"

# recompute summary including scrambled-side mean DMR size
raw <- read.csv("results/dmr_benchmark_readcount_real/readcount_real_per_perm.csv")

real_rows <- raw %>% filter(is_real)
scr_rows  <- raw %>% filter(!is_real)

real_summary <- real_rows %>%
  select(window_size, minDiff, n_real = n_dmrs, size_real_mean = size_mean)

scr_summary <- scr_rows %>%
  group_by(window_size, minDiff) %>%
  summarise(
    n_scr_mean    = mean(n_dmrs, na.rm = TRUE),
    size_scr_mean = mean(size_mean, na.rm = TRUE),  # mean DMR size across the
                                                      # 20 scrambled permutations
    .groups = "drop"
  )

compact <- real_summary %>%
  left_join(scr_summary, by = c("window_size", "minDiff")) %>%
  mutate(
    # coverage = number of DMRs x their actual mean size (Mb), not window
    # size -- the original mistake assumed every DMR is exactly one window
    # wide, which the size_real_mean column shows is false (DMRs merge
    # across adjacent windows, especially at loose thresholds)
    coverage_real_mb = n_real * size_real_mean / 1e6,
    coverage_scr_mb  = n_scr_mean * size_scr_mean / 1e6
  ) %>%
  arrange(window_size, minDiff)

write.csv(compact, file.path(OUT_DIR, "readcount_real_coverage_summary.csv"), row.names = FALSE)
message("Wrote: readcount_real_coverage_summary.csv")
print(compact[, c("window_size","minDiff","n_real","size_real_mean",
                  "n_scr_mean","size_scr_mean","coverage_real_mb","coverage_scr_mb")])

compact$minDiff_lab <- paste0("minDiff = ", compact$minDiff)
compact$minDiff_lab <- factor(compact$minDiff_lab,
  levels = paste0("minDiff = ", c(0.1, 0.2, 0.3, 0.4)))
compact <- compact[!is.na(compact$coverage_real_mb) & !is.na(compact$coverage_scr_mb), ]

cov_long <- rbind(
  data.frame(window_size = compact$window_size * 0.95, mb = compact$coverage_real_mb,
            minDiff_lab = compact$minDiff_lab, series = "Real", stringsAsFactors = FALSE),
  data.frame(window_size = compact$window_size * 1.05, mb = compact$coverage_scr_mb,
            minDiff_lab = compact$minDiff_lab,
            series = "Scrambled (mean of 20 read-count permutations)",
            stringsAsFactors = FALSE)
)
cov_long$series <- factor(cov_long$series,
  levels = c("Real", "Scrambled (mean of 20 read-count permutations)"))

p_cov <- ggplot(cov_long, aes(x = window_size, y = mb, colour = series, shape = series)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.8, fill = "white", stroke = 0.8) +
  scale_colour_manual(values = c("Real" = RED, "Scrambled (mean of 20 read-count permutations)" = GREY)) +
  scale_shape_manual(values = c("Real" = 17, "Scrambled (mean of 20 read-count permutations)" = 2)) +
  scale_x_log10(breaks = c(100, 200, 300, 500, 1000, 2000)) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~minDiff_lab, nrow = 1) +
  labs(
    title = "Read-count permutation: genome coverage (Mb)",
    subtitle = "Coverage = number of DMRs x their actual mean called size, not window size.\nAt loose thresholds, DMRs merge across adjacent windows, so size > window size.",
    x = "Bin/window size (bp)", y = "DMR genome coverage (Mb, log scale)",
    colour = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(OUT_DIR, "readcount_real_coverage.pdf"), p_cov,
       width = 14, height = 5.5, device = cairo_pdf)
message("Saved: readcount_real_coverage.pdf")
