.libPaths("~/R/library")
library(ggplot2)
library(patchwork)

# =============================================================================
# plot_permutations.R  (polished v3 — slide 13 style)
#
# Drop-in replacement. Same input CSV, same output folder, same output filename.
#
# v3 fixes vs v2:
#   - Scrambled now appears in the legend (was missing in v2 output)
#   - Panel B legend has 3 entries: dashed null line, 95% CI, S/N ratio
#     (matches slide 13 exactly)
# =============================================================================

INPUT_CSV <- "results/dmr_benchmark/benchmark_permutations_summary.csv"
PLOT_DIR  <- "results/dmr_benchmark/plots"
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(INPUT_CSV)

# Ratio + 95% CI on the ratio (n = 20 seeds)
n_seeds <- 20
df$se_scrambled  <- df$sd_scrambled / sqrt(n_seeds)
df$ci95_scr      <- qt(0.975, df = n_seeds - 1) * df$se_scrambled
df$ratio         <- df$n_real / df$mean_scrambled
df$ratio_lo      <- df$n_real / (df$mean_scrambled + df$ci95_scr)
df$ratio_hi      <- df$n_real / pmax(df$mean_scrambled - df$ci95_scr, 1e-9)

# ---- Palette ----------------------------------------------------------------
COL_REAL  <- "#2C5F8D"   # navy (real)
COL_SCRAM <- "#B0BEC5"   # light grey-blue (scrambled mean)
COL_RATIO <- "#C0392B"   # brick red (S/N ratio)

# ---- Theme ------------------------------------------------------------------
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
    legend.text        = element_text(size = 10),
    legend.key.width   = unit(1.4, "lines"),
    legend.background  = element_rect(fill = scales::alpha("white", 0.92),
                                      colour = "grey85", linewidth = 0.3),
    legend.margin      = margin(2, 6, 2, 6),
    plot.margin        = margin(10, 14, 10, 14)
  )

x_scale <- scale_x_log10(
  breaks = c(100, 200, 300, 500, 1000, 2000),
  labels = c("100", "200", "300", "500", "1000", "2000"),
  expand = expansion(mult = c(0.05, 0.05))
)

# ---- Panel A: real vs scrambled mean +/- SD ---------------------------------
# Long format + factor on series ensures BOTH appear in the legend
df_A <- rbind(
  data.frame(window_size = df$window_size,
             y           = df$n_real,
             ymin        = df$n_real,
             ymax        = df$n_real,
             series      = "Real DMRs",
             stringsAsFactors = FALSE),
  data.frame(window_size = df$window_size,
             y           = df$mean_scrambled,
             ymin        = df$mean_scrambled - df$sd_scrambled,
             ymax        = df$mean_scrambled + df$sd_scrambled,
             series      = "Scrambled mean \u00b1 SD",
             stringsAsFactors = FALSE)
)
df_A$series <- factor(df_A$series,
                       levels = c("Scrambled mean \u00b1 SD", "Real DMRs"))

pA <- ggplot(df_A, aes(x = window_size, y = y,
                        colour = series, linetype = series, shape = series,
                        group = series)) +
  geom_ribbon(data = df_A[df_A$series == "Scrambled mean \u00b1 SD", ],
              aes(ymin = ymin, ymax = ymax),
              fill = "grey75", alpha = 0.45, inherit.aes = FALSE) +
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
# Three legend entries via dummy columns: null baseline, 95% CI, S/N ratio
df_B <- df
df_B$marker <- "S/N ratio"

pB <- ggplot(df_B, aes(x = window_size)) +
  # 95% CI ribbon (mapped via fill so it gets a legend entry)
  geom_ribbon(aes(ymin = ratio_lo, ymax = ratio_hi, fill = "95% CI"),
              alpha = 0.18) +
  # Null baseline (mapped via linetype = "Null baseline (ratio = 1)")
  geom_hline(aes(yintercept = 1, linetype = "Null baseline (ratio = 1)"),
             colour = "grey60") +
  # S/N ratio line + points
  geom_line(aes(y = ratio, colour = marker), linewidth = 1.0, lineend = "round") +
  geom_point(aes(y = ratio, colour = marker), fill = "white", shape = 21,
             size = 2.8, stroke = 1.0) +
  scale_colour_manual(name = NULL, values = c("S/N ratio" = COL_RATIO)) +
  scale_fill_manual(  name = NULL, values = c("95% CI"    = COL_RATIO)) +
  scale_linetype_manual(name = NULL,
                        values = c("Null baseline (ratio = 1)" = "dashed")) +
  guides(
    linetype = guide_legend(order = 1,
                            override.aes = list(colour = "grey60")),
    fill     = guide_legend(order = 2,
                            override.aes = list(alpha = 0.4)),
    colour   = guide_legend(order = 3,
                            override.aes = list(shape = 21, fill = "white"))
  ) +
  x_scale +
  labs(
    title    = "B   Signal/noise ratio (95% CI)",
    subtitle = "strict threshold, minDiff=0.2",
    x = "Bin size (bp)",
    y = "Signal/noise ratio (real / mean scrambled)"
  ) +
  theme_radu +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.spacing.y = unit(2, "pt"))

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
