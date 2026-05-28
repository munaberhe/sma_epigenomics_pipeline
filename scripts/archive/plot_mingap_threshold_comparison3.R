.libPaths("~/R/library")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
# plot_mingap_threshold_comparison.R  (polished v3 — slide 15 style)
# Drop-in replacement for the original script.
#
# v3 changes:
#   - Stratified scramble dropped (per user)
#   - TWO output PDFs:
#       1. ..._all_series.pdf      \u2014 shows Label-swap and Read-count
#                                    each at minDiff=0.2 and minDiff=0.4
#                                    (4 series total, all visible)
#       2. ..._consolidated.pdf    \u2014 matches slide 15 exactly: groups the
#                                    flat minDiff=0.4 lines into one grey
#                                    "minDiff=0.4 (near-zero)" entry so
#                                    panel A is uncluttered
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

# Drop Stratified
df_all <- df_all %>% filter(scramble_method != "Stratified scramble")
df_all$minDiff <- factor(df_all$minDiff, levels = c("0.2", "0.4 (Radu)"))

# ---- Palette ----------------------------------------------------------------
NAVY       <- "#2C5F8D"
NAVY_LIGHT <- "#7FB0D3"
RED        <- "#C0392B"
RED_LIGHT  <- "#E8A39B"
GREY       <- "#9E9E9E"

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
  breaks = c(100, 200, 500, 1000, 2000),
  labels = c("100", "200", "500", "1000", "2000"),
  expand = expansion(mult = c(0.05, 0.05))
)

# =============================================================================
# OUTPUT 1: All series shown separately
# =============================================================================
df_all$key <- factor(
  paste0(df_all$scramble_method, " \u00b7 minDiff=", df_all$minDiff),
  levels = c(
    "Label-swap \u00b7 minDiff=0.2",             "Label-swap \u00b7 minDiff=0.4 (Radu)",
    "Read count permutation \u00b7 minDiff=0.2", "Read count permutation \u00b7 minDiff=0.4 (Radu)"))

key_cols_all <- setNames(c(NAVY, NAVY_LIGHT, RED, RED_LIGHT), levels(df_all$key))
key_lty_all  <- setNames(c("solid", "longdash", "solid", "longdash"), levels(df_all$key))
key_shp_all  <- setNames(c(15, 0, 18, 5), levels(df_all$key))

df_strict <- df_all %>% filter(mode == "strict", window_size > 0)

pA1 <- ggplot(df_strict,
       aes(x = window_size, y = n_real,
           colour = key, linetype = key, shape = key, group = key)) +
  geom_line(linewidth = 1.0, lineend = "round") +
  geom_point(size = 2.8, stroke = 0.9, fill = "white") +
  scale_colour_manual(values = key_cols_all) +
  scale_linetype_manual(values = key_lty_all) +
  scale_shape_manual(values = key_shp_all) +
  x_scale +
  labs(
    title    = "A   Real DMR counts \u2014 strict threshold",
    subtitle = "minDiff threshold comparison (neighbourhood method)",
    x = "minGap (bp)", y = "n real DMRs"
  ) +
  theme_radu +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1))

df_strict_ratio <- df_strict %>% filter(!is.infinite(ratio), !is.na(ratio))

pB1 <- ggplot(df_strict_ratio,
       aes(x = window_size, y = ratio,
           colour = key, linetype = key, shape = key, group = key)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  annotate("text", x = 100, y = 1, label = "  Null baseline (ratio = 1)",
           hjust = 0, vjust = -0.5, size = 3.3, colour = "grey45",
           fontface = "italic") +
  geom_line(linewidth = 1.0, lineend = "round") +
  geom_point(size = 2.8, stroke = 0.9, fill = "white") +
  scale_colour_manual(values = key_cols_all) +
  scale_linetype_manual(values = key_lty_all) +
  scale_shape_manual(values = key_shp_all) +
  x_scale +
  labs(
    title    = "B   Signal/noise ratio \u2014 strict threshold",
    subtitle = "minDiff threshold comparison",
    x = "minGap (bp)", y = "Signal/noise ratio (real / scrambled)"
  ) +
  theme_radu +
  theme(legend.position = "none")

combined_all <- (pA1 | pB1) +
  plot_annotation(
    caption = "Solid lines = minDiff=0.2; dashed lines = minDiff=0.4. minDiff=0.4 recovers near-zero DMRs under strict threshold across all window sizes and null models.",
    theme = theme(plot.caption = element_text(face = "italic",
                                              colour = "grey40", hjust = 0))
  )

ggsave(file.path(PLOT_DIR, "nb_mingap_threshold_comparison_all_series.pdf"),
       combined_all, width = 14, height = 6,
       device = grDevices::cairo_pdf)
message("Saved: nb_mingap_threshold_comparison_all_series.pdf")

# =============================================================================
# OUTPUT 2: Consolidated (slide 15 style — groups flat 0.4 lines into one)
# =============================================================================
# Build a consolidated key: 0.2 series remain as-is; both 0.4 series get
# merged into "minDiff=0.4 (near-zero)" plotted as a single grey series
# (we average them so the line position is meaningful at zero anyway).
df_cons <- df_strict %>%
  mutate(key_cons = case_when(
    minDiff == "0.2"        ~ paste0("minDiff=0.2 \u00b7 ", scramble_method),
    minDiff == "0.4 (Radu)" ~ "minDiff=0.4 (near-zero)",
    TRUE                    ~ as.character(scramble_method)
  ))

# For panel A, average the 0.4 lines into one (they're both flat at zero)
df_cons_A <- df_cons %>%
  group_by(window_size, key_cons) %>%
  summarise(n_real = mean(n_real, na.rm = TRUE), .groups = "drop")

# For panel B, keep both 0.4 lines (they DO differ on the ratio scale)
df_cons_B <- df_strict_ratio

key_cons_levels <- c(
  "minDiff=0.2 \u00b7 Label-swap",
  "minDiff=0.2 \u00b7 Read count permutation",
  "minDiff=0.4 (near-zero)"
)
df_cons_A$key_cons <- factor(df_cons_A$key_cons, levels = key_cons_levels)

cons_cols <- setNames(c(NAVY, RED, GREY), key_cons_levels)
cons_lty  <- setNames(c("solid", "solid", "dashed"), key_cons_levels)
cons_shp  <- setNames(c(15, 18, 22), key_cons_levels)

pA2 <- ggplot(df_cons_A,
       aes(x = window_size, y = n_real,
           colour = key_cons, linetype = key_cons, shape = key_cons,
           group = key_cons)) +
  geom_line(linewidth = 1.0, lineend = "round") +
  geom_point(size = 2.8, stroke = 0.9, fill = "white") +
  scale_colour_manual(values = cons_cols) +
  scale_linetype_manual(values = cons_lty) +
  scale_shape_manual(values = cons_shp) +
  x_scale +
  labs(
    title    = "A   Real DMR counts \u2014 strict threshold",
    subtitle = "minDiff threshold comparison (neighbourhood method)",
    x = "minGap (bp)", y = "n real DMRs"
  ) +
  theme_radu +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1))

# Panel B for consolidated: keep all 4 series so ratio differences are visible
df_cons_B$key_cons <- df_cons_B$key
key_cols_B <- key_cols_all
key_lty_B  <- key_lty_all
key_shp_B  <- key_shp_all

pB2 <- ggplot(df_cons_B,
       aes(x = window_size, y = ratio,
           colour = key_cons, linetype = key_cons, shape = key_cons,
           group = key_cons)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  annotate("text", x = 100, y = 1, label = "  Null baseline (ratio = 1)",
           hjust = 0, vjust = -0.5, size = 3.3, colour = "grey45",
           fontface = "italic") +
  geom_line(linewidth = 1.0, lineend = "round") +
  geom_point(size = 2.8, stroke = 0.9, fill = "white") +
  scale_colour_manual(values = key_cols_B) +
  scale_linetype_manual(values = key_lty_B) +
  scale_shape_manual(values = key_shp_B) +
  x_scale +
  labs(
    title    = "B   Signal/noise ratio \u2014 strict threshold",
    subtitle = "minDiff threshold comparison",
    x = "minGap (bp)", y = "Signal/noise ratio (real / scrambled)"
  ) +
  theme_radu +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1))

combined_cons <- (pA2 | pB2) +
  plot_annotation(
    caption = "Solid lines = minDiff=0.2; dashed lines = minDiff=0.4. minDiff=0.4 recovers near-zero DMRs under strict threshold across all window sizes and null models.",
    theme = theme(plot.caption = element_text(face = "italic",
                                              colour = "grey40", hjust = 0))
  )

ggsave(file.path(PLOT_DIR, "nb_mingap_threshold_comparison_consolidated.pdf"),
       combined_cons, width = 14, height = 6,
       device = grDevices::cairo_pdf)
message("Saved: nb_mingap_threshold_comparison_consolidated.pdf")
