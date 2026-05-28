.libPaths("~/R/library")
library(ggplot2)
library(patchwork)

# =============================================================================
# plot_permutations.R  (polished — style-matched to slide 13)
#
# Drop-in replacement. Same input CSV, same output folder, same output PDF
# filename. Only the plotting layer is rewritten.
#
# Input:
#   results/dmr_benchmark/benchmark_permutations_summary.csv
#     columns: window_size, n_real, mean_scrambled, sd_scrambled, ...
#
# Output:
#   results/dmr_benchmark/plots/radu_panels_read_count_permutation_20seeds.pdf
#
# Key style fixes vs the old version:
#   1. Bold A/B panel-title prefix; italic grey caption footer
#   2. Real line in navy (solid, square markers), scrambled mean in light
#      navy (dashed, open square markers) so they read as the SAME method
#      but distinct series — matches slide 13
#   3. Panel B adds a signal/noise ratio panel with 95% CI shading (this is
#      what slide 13's panel B actually shows — your old script plotted
#      delta coverage instead)
#   4. Log-scaled x-axis so 100/200/300 don't squash
#   5. "Null baseline (ratio = 1)" annotation on dashed reference line
# =============================================================================

INPUT_CSV <- "results/dmr_benchmark/benchmark_permutations_summary.csv"
PLOT_DIR  <- "results/dmr_benchmark/plots"
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(INPUT_CSV)

# Coverage in Mb (kept for caption-only reference; plot uses raw counts)
df$cov_real_mb     <- df$n_real         * df$window_size / 1e6
df$cov_scr_mean_mb <- df$mean_scrambled * df$window_size / 1e6
df$cov_scr_sd_mb   <- df$sd_scrambled   * df$window_size / 1e6
df$cov_diff_mb     <- df$cov_real_mb - df$cov_scr_mean_mb

# Ratio + 95% CI on the ratio (using SE = SD/sqrt(n), n = 20 seeds)
n_seeds <- 20
df$se_scrambled  <- df$sd_scrambled / sqrt(n_seeds)
df$ci95_scr      <- qt(0.975, df = n_seeds - 1) * df$se_scrambled
df$ratio         <- df$n_real / df$mean_scrambled
df$ratio_lower      <- df$n_real / (df$mean_scrambled + df$ci95_scr)
df$ratio_upper      <- df$n_real / pmax(df$mean_scrambled - df$ci95_scr, 1e-9)

# ---- Palette / theme --------------------------------------------------------
COL_REAL  <- "#2C5F8D"   # navy (solid)
COL_SCRAM <- "#7FB0D3"   # light navy (dashed)
COL_RATIO <- "#C0392B"   # brick red

theme_radu <- theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 13, hjust = 0,
                                      margin = margin(b = 4)),
    plot.subtitle      = element_text(size = 12, hjust = 0,
                                      margin = margin(b = 8)),
    plot.caption       = element_text(face = "italic", colour = "grey40",
                                      hjust = 0, size = 10,
                                      margin = margin(t = 10)),
    plot.caption.position = "plot",
    axis.title         = element_text(size = 12),
    axis.text          = element_text(size = 11, colour = "grey20"),
    axis.line          = element_line(colour = "grey30", linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.major.x = element_line(colour = "grey95", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11),
    legend.key.width   = unit(1.4, "lines"),
    legend.background  = element_rect(fill = scales::alpha("white", 0.9),
                                      colour = "grey85", linewidth = 0.3),
    legend.margin      = margin(2, 6, 2, 6),
    plot.margin        = margin(10, 14, 10, 14)
  )

x_scale <- scale_x_log10(
  breaks = c(100, 200, 300, 500, 1000, 2000),
  labels = c("100", "200", "300", "500", "1000", "2000"),
  expand = expansion(mult = c(0.05, 0.05))
)

# ---- Panel A: DMR counts, real vs scrambled mean +/- SD --------------------
# Reshape so both series go through scale_colour_manual + scale_linetype_manual
# in one mapping — this is what makes BOTH appear in the legend.
df_long <- rbind(
  data.frame(window_size = df$window_size, y = df$n_real,
             ymin = df$n_real, ymax = df$n_real,
             series = "Real DMRs", stringsAsFactors = FALSE),
  data.frame(window_size = df$window_size, y = df$mean_scrambled,
             ymin = df$mean_scrambled - df$sd_scrambled,
             ymax = df$mean_scrambled + df$sd_scrambled,
             series = "Scrambled mean \u00b1 SD", stringsAsFactors = FALSE)
)
df_long$series <- factor(df_long$series,
                          levels = c("Real DMRs", "Scrambled mean \u00b1 SD"))

pA <- ggplot(df_long, aes(x = window_size, y = y,
                           colour = series, linetype = series, shape = series,
                           group = series)) +
  # SD ribbon under the scrambled line only (df_long carries y/ymin/ymax)
  geom_ribbon(data = df_long[df_long$series == "Scrambled mean \u00b1 SD", ],
              aes(x = window_size, ymin = ymin, ymax = ymax),
              fill = "grey75", alpha = 0.45,
              inherit.aes = FALSE) +
  geom_line(linewidth = 1.0, lineend = "round") +
  geom_point(size = 2.8, stroke = 0.9, fill = "white") +
  scale_colour_manual(values = c("Real DMRs"               = COL_REAL,
                                  "Scrambled mean \u00b1 SD" = COL_SCRAM)) +
  scale_linetype_manual(values = c("Real DMRs"               = "solid",
                                    "Scrambled mean \u00b1 SD" = "longdash")) +
  scale_shape_manual(values = c("Real DMRs"               = 15,
                                 "Scrambled mean \u00b1 SD" = 0)) +
  x_scale +
  labs(
    title    = "A   DMR counts \u2014 real vs null (20 permutations)",
    subtitle = "null model: Read count permutation, DMRcaller-B (strict)",
    x = "Bin size (bp)",
    y = "Number of DMRs"
  ) +
  theme_radu +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1))

# ---- Panel B: signal/noise ratio with 95% CI -------------------------------
pB <- ggplot(df, aes(x = window_size, y = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(x = window_size, ymin = ratio_lower, ymax = ratio_upper),
              fill = COL_RATIO, alpha = 0.18) +
  geom_line(colour = COL_RATIO, linewidth = 1.0) +
  geom_point(colour = COL_RATIO, fill = "white", shape = 21,
             size = 2.8, stroke = 1.0) +
  annotate("text", x = 100, y = 1, label = "  Null baseline (ratio = 1)",
           hjust = 0, vjust = -0.5, size = 3.3, colour = "grey45",
           fontface = "italic") +
  x_scale +
  labs(
    title    = "B   Signal/noise ratio (95% CI)",
    subtitle = "strict threshold, minDiff=0.2",
    x = "Bin size (bp)",
    y = "Signal/noise ratio (real / mean scrambled)"
  ) +
  theme_radu +
  theme(legend.position = "none")

combined <- (pA | pB) +
  plot_annotation(
    caption = "Shaded band = mean \u00b1 1 SD (panel A) and 95% CI (panel B) across 20 permutation seeds.",
    theme = theme(plot.caption = element_text(face = "italic",
                                              colour = "grey40", hjust = 0))
  )

outfile <- file.path(PLOT_DIR, "radu_panels_read_count_permutation_20seeds.pdf")
ggsave(outfile, combined, width = 14, height = 5.5,
       device = grDevices::cairo_pdf)
message("Saved: ", basename(outfile))
