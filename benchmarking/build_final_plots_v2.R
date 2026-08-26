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
SELECTED_WS <- 300    # selected pipeline bin size
SELECTED_MD <- 0.2    # selected pipeline minProportionDifference

# PLOT 1: label-swap, real replicate-partition permutation (wide format
# for ribbon shading)
ls_df <- read.csv("results/dmr_benchmark_labelswap_real/labelswap_real_summary.csv")
ls_df <- ls_df[ls_df$method == "bins", ]
ls_df <- ls_df[!is.na(ls_df$n_real), ]
ls_df$window_size <- as.numeric(ls_df$window_size)
ls_df <- ls_df[order(ls_df$window_size), ]

p_ls <- ggplot(ls_df, aes(x = window_size)) +
  geom_ribbon(aes(ymin = pmin(n_real, n_scr_mean), ymax = pmax(n_real, n_scr_mean)),
              fill = NAVY, alpha = 0.12) +
  geom_line(aes(y = n_real, colour = "Real"), linewidth = 1.0) +
  geom_point(aes(y = n_real, colour = "Real", shape = "Real"), size = 3, fill = "white", stroke = 0.9) +
  geom_line(aes(y = n_scr_mean, colour = "Scrambled (mean of 9 real partitions)"), linewidth = 1.0) +
  geom_point(aes(y = n_scr_mean, colour = "Scrambled (mean of 9 real partitions)",
                shape = "Scrambled (mean of 9 real partitions)"), size = 3, fill = "white", stroke = 0.9) +
  geom_vline(xintercept = SELECTED_WS, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  annotate("text", x = SELECTED_WS, y = max(ls_df$n_real, na.rm = TRUE) * 1.5,
           label = paste0("Selected bin size (", SELECTED_WS, "bp)"),
           angle = 90, vjust = -0.5, size = 3.2, colour = "grey30") +
  scale_colour_manual(values = c("Real" = NAVY, "Scrambled (mean of 9 real partitions)" = GREY),
                      name = NULL) +
  scale_shape_manual(values = c("Real" = 15, "Scrambled (mean of 9 real partitions)" = 0),
                     name = NULL) +
  scale_x_log10(breaks = c(100, 200, 300, 500, 1000, 2000)) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Label-swap: real replicate-partition permutation",
    subtitle = "DMRcaller-B, strict thresholds. Null = all 9 non-identity 3-vs-3 replicate\npartitions of the 6 ASO_VPA/ASO_CTRL replicates (exhaustive, not sampled).\nDashed line marks the selected pipeline bin size.",
    x = "Bin/window size (bp)", y = "Number of DMRs (log scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "labelswap_real_corrected.pdf"), p_ls,
       width = 9, height = 6.5, device = cairo_pdf)
message("Saved: labelswap_real_corrected.pdf")

# PLOT 2: read-count, DMR count, faceted by threshold, selected panel highlighted
rc_df <- read.csv("results/dmr_benchmark_readcount_real/readcount_real_summary.csv")
rc_df <- rc_df[!is.na(rc_df$n_real) & !is.na(rc_df$n_scr_mean), ]
rc_df$window_size <- as.numeric(rc_df$window_size)
rc_df <- rc_df[order(rc_df$window_size), ]
rc_df$is_selected <- rc_df$minDiff == SELECTED_MD
rc_df$minDiff_lab <- paste0("minDiff = ", rc_df$minDiff,
                            ifelse(rc_df$is_selected, "  (selected)", ""))
rc_df$minDiff_lab <- factor(rc_df$minDiff_lab,
  levels = unique(rc_df$minDiff_lab[order(rc_df$minDiff)]))

p_rc <- ggplot(rc_df, aes(x = window_size)) +
  geom_rect(data = rc_df[rc_df$is_selected, ][1, ],
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "#FFF6E5", inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = pmin(n_real, n_scr_mean), ymax = pmax(n_real, n_scr_mean)),
              fill = RED, alpha = 0.10) +
  geom_line(aes(y = n_real, colour = "Real"), linewidth = 1.0) +
  geom_point(aes(y = n_real, colour = "Real", shape = "Real"), size = 2.6, fill = "white", stroke = 0.8) +
  geom_line(aes(y = n_scr_mean, colour = "Scrambled (mean of 20 read-count permutations)"), linewidth = 1.0) +
  geom_point(aes(y = n_scr_mean, colour = "Scrambled (mean of 20 read-count permutations)",
                shape = "Scrambled (mean of 20 read-count permutations)"), size = 2.6, fill = "white", stroke = 0.8) +
  geom_vline(xintercept = SELECTED_WS, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = c("Real" = RED, "Scrambled (mean of 20 read-count permutations)" = GREY),
                      name = NULL) +
  scale_shape_manual(values = c("Real" = 17, "Scrambled (mean of 20 read-count permutations)" = 2),
                     name = NULL) +
  scale_x_log10(breaks = c(100, 200, 300, 500, 1000, 2000)) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~minDiff_lab, nrow = 1) +
  labs(
    title = "Read-count permutation: threshold sweep (DMR count)",
    subtitle = "Shaded panel = selected threshold (minDiff = 0.2). Dashed line = selected bin size (300bp).",
    x = "Bin/window size (bp)", y = "Number of DMRs (log scale)"
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

# PLOT 3: read-count, genome coverage (Mb), same shading/highlighting
raw <- read.csv("results/dmr_benchmark_readcount_real/readcount_real_per_perm.csv")
real_rows <- raw %>% filter(is_real)
scr_rows  <- raw %>% filter(!is_real)

real_summary <- real_rows %>%
  select(window_size, minDiff, n_real = n_dmrs, size_real_mean = size_mean)
scr_summary <- scr_rows %>%
  group_by(window_size, minDiff) %>%
  summarise(n_scr_mean = mean(n_dmrs, na.rm = TRUE),
            size_scr_mean = mean(size_mean, na.rm = TRUE), .groups = "drop")

cov_df <- real_summary %>%
  left_join(scr_summary, by = c("window_size", "minDiff")) %>%
  mutate(coverage_real_mb = n_real * size_real_mean / 1e6,
         coverage_scr_mb  = n_scr_mean * size_scr_mean / 1e6) %>%
  filter(!is.na(coverage_real_mb), !is.na(coverage_scr_mb)) %>%
  arrange(window_size, minDiff)

cov_df$window_size <- as.numeric(cov_df$window_size)
cov_df$is_selected <- cov_df$minDiff == SELECTED_MD
cov_df$minDiff_lab <- paste0("minDiff = ", cov_df$minDiff,
                             ifelse(cov_df$is_selected, "  (selected)", ""))
cov_df$minDiff_lab <- factor(cov_df$minDiff_lab,
  levels = unique(cov_df$minDiff_lab[order(cov_df$minDiff)]))

p_cov <- ggplot(cov_df, aes(x = window_size)) +
  geom_rect(data = cov_df[cov_df$is_selected, ][1, ],
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "#FFF6E5", inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = pmin(coverage_real_mb, coverage_scr_mb),
                  ymax = pmax(coverage_real_mb, coverage_scr_mb)),
              fill = RED, alpha = 0.10) +
  geom_line(aes(y = coverage_real_mb, colour = "Real"), linewidth = 1.0) +
  geom_point(aes(y = coverage_real_mb, colour = "Real", shape = "Real"), size = 2.6, fill = "white", stroke = 0.8) +
  geom_line(aes(y = coverage_scr_mb, colour = "Scrambled (mean of 20 read-count permutations)"), linewidth = 1.0) +
  geom_point(aes(y = coverage_scr_mb, colour = "Scrambled (mean of 20 read-count permutations)",
                shape = "Scrambled (mean of 20 read-count permutations)"), size = 2.6, fill = "white", stroke = 0.8) +
  geom_vline(xintercept = SELECTED_WS, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = c("Real" = RED, "Scrambled (mean of 20 read-count permutations)" = GREY),
                      name = NULL) +
  scale_shape_manual(values = c("Real" = 17, "Scrambled (mean of 20 read-count permutations)" = 2),
                     name = NULL) +
  scale_x_log10(breaks = c(100, 200, 300, 500, 1000, 2000)) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~minDiff_lab, nrow = 1) +
  labs(
    title = "Read-count permutation: genome coverage (Mb)",
    subtitle = "Coverage = DMR count x actual mean called size. Shaded panel = selected threshold.\nDashed line = selected bin size (300bp).",
    x = "Bin/window size (bp)", y = "DMR genome coverage (Mb, log scale)"
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

message("\nAll three plots saved in: ", OUT_DIR)
