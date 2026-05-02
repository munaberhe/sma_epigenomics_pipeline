.libPaths("~/R/library")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
# plot_mingap_threshold_comparison.R  (polished — style-matched to slide 15)
#
# Drop-in replacement. Same inputs (CSV + RDS checkpoint), same output folder,
# same output PDF filename.
#
# Key style fixes vs the old version:
#   1. Bold A/B panel-title prefix; italic grey caption footer
#   2. Solid lines = minDiff=0.2; dashed lines = minDiff=0.4 (Radu's params)
#      drawn on the SAME panel rather than facetted, matching slide 15
#   3. Navy / red / green palette for the three null models
#   4. Filled markers (square/diamond/triangle) with white centres
#   5. Log-scaled x-axis so 100/200 don't squash
#   6. "Null baseline (ratio = 1)" annotation on dashed reference line
# =============================================================================

OUT_DIR  <- "results/dmr_benchmark"
PLOT_DIR <- "results/dmr_benchmark/plots_v2"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- minDiff=0.2 ------------------------------------------------------------
ls <- read.csv(file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap.csv"))
ls <- ls %>% filter(window_size > 0) %>%
  mutate(scramble_method = "label_swap", minDiff = "0.2")

ckpt <- readRDS(file.path(OUT_DIR, "checkpoint_neighbourhood_mingap_v2.rds"))
strat <- do.call(rbind, ckpt) %>%
  rename(window_size = window_size) %>%
  mutate(minDiff = "0.2")

df_02 <- bind_rows(
  ls    %>% select(window_size, mode, n_real, n_scrambled, ratio, scramble_method, minDiff),
  strat %>% select(window_size, mode, n_real, n_scrambled, ratio, scramble_method, minDiff)
)

# ---- minDiff=0.4 (Radu) -----------------------------------------------------
fin <- read.csv(file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap_final.csv"))
df_04 <- fin %>%
  filter(minsize == 100) %>%
  rename(window_size = mingap) %>%
  mutate(minDiff = "0.4 (Radu)", mode = "strict") %>%
  select(window_size, mode, n_real, n_scrambled, ratio, scramble_method, minDiff)

df_all <- bind_rows(df_02, df_04)

df_all$scramble_method <- recode(df_all$scramble_method,
  label_swap = "Label-swap",
  stratified = "Stratified scramble",
  archie     = "Read count permutation")

df_all$minDiff <- factor(df_all$minDiff, levels = c("0.2", "0.4 (Radu)"))

# Build a combined "key" combining null model + minDiff so we get one legend
# entry per series and can vary linetype within colour
df_all$key <- factor(
  paste0(df_all$scramble_method, " \u00b7 minDiff=", df_all$minDiff),
  levels = c(
    "Label-swap \u00b7 minDiff=0.2",             "Label-swap \u00b7 minDiff=0.4 (Radu)",
    "Stratified scramble \u00b7 minDiff=0.2",    "Stratified scramble \u00b7 minDiff=0.4 (Radu)",
    "Read count permutation \u00b7 minDiff=0.2", "Read count permutation \u00b7 minDiff=0.4 (Radu)"))

# ---- Palette ----------------------------------------------------------------
col_null <- c(
  "Label-swap"             = "#2C5F8D",   # navy
  "Stratified scramble"    = "#2E7D32",   # green
  "Read count permutation" = "#C0392B")   # brick red
col_null_light <- c(
  "Label-swap"             = "#7FB0D3",
  "Stratified scramble"    = "#9CCBA0",
  "Read count permutation" = "#E8A39B")

key_cols <- setNames(
  c(col_null["Label-swap"],             col_null_light["Label-swap"],
    col_null["Stratified scramble"],    col_null_light["Stratified scramble"],
    col_null["Read count permutation"], col_null_light["Read count permutation"]),
  levels(df_all$key))

key_lty <- setNames(
  c("solid", "dashed", "solid", "dashed", "solid", "dashed"),
  levels(df_all$key))

key_shp <- setNames(
  c(15, 0,  17, 2,  18, 5),
  levels(df_all$key))

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
    legend.background  = element_rect(fill = scales::alpha("white", 0.9),
                                      colour = "grey85", linewidth = 0.3),
    legend.margin      = margin(2, 6, 2, 6),
    plot.margin        = margin(10, 14, 10, 14)
  )

x_scale <- scale_x_log10(
  breaks = c(100, 200, 500, 1000, 2000),
  labels = c("100", "200", "500", "1000", "2000"),
  expand = expansion(mult = c(0.05, 0.05))
)

# ---- Panel A: real DMR counts ----------------------------------------------
df_strict <- df_all %>% filter(mode == "strict", window_size > 0)

pA <- ggplot(df_strict,
       aes(x = window_size, y = n_real,
           colour = key, linetype = key, shape = key, group = key)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8, stroke = 0.8, fill = "white") +
  scale_colour_manual(values = key_cols) +
  scale_linetype_manual(values = key_lty) +
  scale_shape_manual(values = key_shp) +
  x_scale +
  labs(
    title    = "A   Real DMR counts \u2014 strict threshold",
    subtitle = "minDiff threshold comparison (neighbourhood method)",
    x = "minGap (bp)",
    y = "n real DMRs"
  ) +
  theme_radu +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1))

# ---- Panel B: signal/noise ratio -------------------------------------------
df_strict_ratio <- df_strict %>%
  filter(!is.infinite(ratio), !is.na(ratio))

pB <- ggplot(df_strict_ratio,
       aes(x = window_size, y = ratio,
           colour = key, linetype = key, shape = key, group = key)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  annotate("text", x = 100, y = 1, label = "  Null baseline (ratio = 1)",
           hjust = 0, vjust = -0.5, size = 3.3, colour = "grey45",
           fontface = "italic") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8, stroke = 0.8, fill = "white") +
  scale_colour_manual(values = key_cols) +
  scale_linetype_manual(values = key_lty) +
  scale_shape_manual(values = key_shp) +
  x_scale +
  labs(
    title    = "B   Signal/noise ratio \u2014 strict threshold",
    subtitle = "minDiff threshold comparison",
    x = "minGap (bp)",
    y = "Signal/noise ratio (real / scrambled)"
  ) +
  theme_radu +
  theme(legend.position = "none")

combined <- (pA | pB) +
  plot_annotation(
    caption = "Solid lines = minDiff=0.2; dashed lines = minDiff=0.4. minDiff=0.4 recovers near-zero DMRs under strict threshold across all window sizes and null models.",
    theme = theme(plot.caption = element_text(face = "italic",
                                              colour = "grey40", hjust = 0))
  )

ggsave(file.path(PLOT_DIR, "nb_mingap_threshold_comparison.pdf"),
       combined, width = 14, height = 6,
       device = grDevices::cairo_pdf)
message("Saved: nb_mingap_threshold_comparison.pdf")
