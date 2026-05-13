.libPaths("~/R/library")
# plot_mingap_results.R  (polished — style-matched to slide 14)
# SMA Epigenomics Pipeline — Muna Berhe, QMUL
#
# Drop-in replacement. Same input CSV, same output folder, same output PDF
# filename. Only the plotting layer is rewritten.
#
# Input:
#   results/dmr_benchmark/parameter_benchmark_neighbourhood_mingap_final.csv
#     columns: mingap, minsize, scramble_method, n_real, n_scrambled, ...
#
# Output:
#   results/dmr_benchmark/plots_v2/mingap_benchmark_counts.pdf
#
# Key style fixes vs the old version:
#   1. Bold A/B panel-title prefix; italic grey caption footer
#   2. Real solid + Scrambled dashed in the SAME panel (matches slide 14
#      which puts both series on one chart, not separate panels)
#   3. Scrambled lines use a lighter shade so they remain visible when
#      sitting on top of the real line
#   4. Navy/red palette matching the rest of the deck
#   5. Filled markers (square/triangle/diamond) with white centres
#   6. Log-scaled x-axis

library(ggplot2)
library(patchwork)
library(scales)

OUT_DIR  <- "results/dmr_benchmark"
PLOT_DIR <- "results/dmr_benchmark/plots_v2"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap_final.csv"))

# Keep mingap numeric for log axis (the old script factored it which forced
# linear spacing and made the small bins squash together)
df$mingap_num <- as.numeric(as.character(df$mingap))

df$minsize <- factor(df$minsize, labels = c("minSize=100", "minSize=200"))
df$scramble_method <- factor(df$scramble_method,
  levels = c("label_swap", "archie", "stratified"),
  labels = c("Label-swap", "Read count permutation", "Stratified scramble"))

# ---- Palette ----------------------------------------------------------------
# Bright = Real, light = Scrambled. Same hue per null model so they read as
# "this null, real" vs "this null, scrambled".
col_real <- c(
  "Label-swap"             = "#2C5F8D",   # navy
  "Read count permutation" = "#C0392B",   # brick red
  "Stratified scramble"    = "#2E7D32")   # green
col_scram <- c(
  "Label-swap"             = "#7FB0D3",
  "Read count permutation" = "#E8A39B",
  "Stratified scramble"    = "#9CCBA0")

shapes_null <- c(
  "Label-swap"             = 15,   # filled square
  "Read count permutation" = 18,   # filled diamond
  "Stratified scramble"    = 17)   # filled triangle

# Long format: one row per (mingap, minsize, null, series)
df_long <- rbind(
  data.frame(mingap_num = df$mingap_num,
             minsize    = df$minsize,
             null       = df$scramble_method,
             n          = df$n_real,
             series     = "real",
             stringsAsFactors = FALSE),
  data.frame(mingap_num = df$mingap_num,
             minsize    = df$minsize,
             null       = df$scramble_method,
             n          = df$n_scrambled,
             series     = "scrambled",
             stringsAsFactors = FALSE)
)

df_long$null <- factor(df_long$null,
  levels = c("Label-swap", "Read count permutation", "Stratified scramble"))

# Build a combined "key" for legend ordering: null × series
df_long$key <- factor(
  paste0(df_long$null, " (", df_long$series, ")"),
  levels = c(
    "Label-swap (real)",             "Label-swap (scrambled)",
    "Read count permutation (real)", "Read count permutation (scrambled)",
    "Stratified scramble (real)",    "Stratified scramble (scrambled)"))

key_cols <- setNames(
  c(col_real["Label-swap"],             col_scram["Label-swap"],
    col_real["Read count permutation"], col_scram["Read count permutation"],
    col_real["Stratified scramble"],    col_scram["Stratified scramble"]),
  levels(df_long$key))
key_lty <- setNames(
  c("solid","dashed","solid","dashed","solid","dashed"),
  levels(df_long$key))
key_shp <- setNames(
  c(15, 0,  18, 5,  17, 2),
  levels(df_long$key))

# Slight x offset so scrambled dashes peek out from under the real line
df_long$x_plot <- ifelse(df_long$series == "real",
                          df_long$mingap_num * 0.975,
                          df_long$mingap_num * 1.025)

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
    strip.background   = element_rect(fill = "grey95", colour = NA),
    strip.text         = element_text(face = "bold", size = 11,
                                      margin = margin(4, 4, 4, 4)),
    legend.title       = element_blank(),
    legend.text        = element_text(size = 10),
    legend.key.width   = unit(1.4, "lines"),
    legend.background  = element_rect(fill = scales::alpha("white", 0.9),
                                      colour = "grey85", linewidth = 0.3),
    legend.margin      = margin(2, 6, 2, 6),
    plot.margin        = margin(10, 14, 10, 14)
  )

x_scale <- scale_x_log10(
  breaks = sort(unique(df$mingap_num)),
  labels = sort(unique(df$mingap_num)),
  expand = expansion(mult = c(0.05, 0.05))
)

# ---- Combined Real-vs-Scrambled plot ---------------------------------------
p <- ggplot(df_long,
       aes(x = x_plot, y = n,
           colour = key, linetype = key, shape = key,
           group = interaction(key, minsize))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6, stroke = 0.8, fill = "white") +
  facet_wrap(~minsize, ncol = 2) +
  scale_colour_manual(values = key_cols) +
  scale_linetype_manual(values = key_lty) +
  scale_shape_manual(values = key_shp) +
  x_scale +
  labs(
    title    = "Neighbourhood method \u2014 Real vs Scrambled DMR Counts",
    subtitle = "strict threshold, minDiff=0.2",
    x = "minGap (bp)",
    y = "Number of DMRs",
    caption = "DMRcaller neighbourhood method. Strict threshold, minDiff=0.2. Scrambled lines offset on x-axis for visibility."
  ) +
  theme_radu +
  theme(legend.position = "top",
        legend.box = "horizontal")

ggsave(file.path(PLOT_DIR, "mingap_benchmark_counts.pdf"),
       p, width = 12, height = 6,
       device = grDevices::cairo_pdf)
message("Saved: mingap_benchmark_counts.pdf")
