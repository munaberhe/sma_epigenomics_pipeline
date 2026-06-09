.libPaths("~/R/library")
library(ggplot2)
library(patchwork)

# =============================================================================
# plot_radu_panels.R  (polished v3 — slide 11/12 style, unified legend)
#
# Drop-in replacement. Same input CSVs, same OUT_DIR/PLOT_DIR, same output
# filenames. Only the plotting layer changed.
#
# v3 fixes vs v2:
#   - Single unified legend per figure (was: two stacked "Method" legends on
#     the right of the patchwork output)
#   - 4 explicit legend entries: B (real), B (scrambled), NF (real), NF (scrambled)
#   - 5% x-offset + longdash so scrambled is clearly separated from real
#   - DMRcaller-NB excluded from BOTH coverage panels AND the panel B legend
# =============================================================================

USE_TURBO <- FALSE
OUT_DIR   <- if (USE_TURBO) "results/dmr_benchmark_turbo" else "results/dmr_benchmark"
PLOT_DIR  <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
message("Reading CSVs from: ", OUT_DIR)

csv_files <- list.files(OUT_DIR, pattern = "parameter_benchmark.*\\.csv",
                        full.names = TRUE)
if (length(csv_files) == 0) stop("No parameter_benchmark*.csv found in ", OUT_DIR)
csv_files <- csv_files[grepl("label_swap", csv_files)]
if (length(csv_files) == 0) stop("No label_swap CSV found in ", OUT_DIR)
df <- do.call(rbind, lapply(csv_files, read.csv))
df$ratio[is.infinite(df$ratio)] <- NA
df$window_size_num <- as.numeric(as.character(df$window_size))

df$method <- factor(df$method,
  levels = c("bins", "neighbourhood", "noise_filter"),
  labels = c("DMRcaller-B", "DMRcaller-NB", "DMRcaller-NF"))

df$scramble_method <- factor(df$scramble_method,
  levels = c("label_swap", "archie_scramble", "stratified_scramble"),
  labels = c("Label-swap", "Read count permutation", "Stratified scramble"))

# Strict mode + Gaussian-only NF
df_s <- df[df$mode == "strict", ]
df_s <- df_s[!(df_s$method == "DMRcaller-NF" &
               df_s$kernel %in% c("uniform", "epanechnicov")), ]

df_s$cov_real_mb      <- df_s$n_real      * df_s$window_size_num / 1e6
df_s$cov_scrambled_mb <- df_s$n_scrambled * df_s$window_size_num / 1e6
df_s$cov_diff_mb      <- df_s$cov_real_mb - df_s$cov_scrambled_mb

# Optional contrast support (no-op for current data which has only ASO)
HAS_CONTRAST <- "contrast" %in% colnames(df_s)
if (HAS_CONTRAST) {
  contrast_levels <- unique(as.character(df_s$contrast))
  contrast_labels <- c(
    VPA_effect = "Scramble_VPA vs Scramble_CTRL",
    ASO_effect = "ASO_VPA vs ASO_CTRL",
    Combined   = "ASO_VPA vs Scramble_CTRL"
  )
} else {
  df_s$contrast   <- "ASO_effect"
  contrast_levels <- "ASO_effect"
  contrast_labels <- c(ASO_effect = "ASO_VPA vs ASO_CTRL")
}

# ---- Palette (slide 11/12) --------------------------------------------------
NAVY        <- "#2C5F8D"
NAVY_LIGHT  <- "#7FB0D3"
RED         <- "#C0392B"
RED_LIGHT   <- "#E8A39B"

# Legend entries: 4 keys (B real/scrambled, NF real/scrambled)
key_levels <- c("DMRcaller-B (real)",  "DMRcaller-B (scrambled)",
                "DMRcaller-NF (real)", "DMRcaller-NF (scrambled)")
key_cols <- setNames(c(NAVY, NAVY_LIGHT, RED, RED_LIGHT), key_levels)
key_lty  <- setNames(c("solid", "longdash", "solid", "longdash"), key_levels)
key_shp  <- setNames(c(15, 0, 18, 5), key_levels)

# Panel B legend (delta — only real method colour, not series)
delta_cols <- c("DMRcaller-B" = NAVY, "DMRcaller-NF" = RED)
delta_shp  <- c("DMRcaller-B" = 15,   "DMRcaller-NF" = 18)

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

# ---- Plotting function ------------------------------------------------------
make_panels <- function(null_label, filename_tag, contrast_key) {
  d <- df_s[df_s$scramble_method == null_label &
            df_s$contrast        == contrast_key &
            df_s$method          != "DMRcaller-NB", ]   # NB excluded
  if (nrow(d) == 0) {
    message("  [skip] no rows for ", null_label, " / ", contrast_key)
    return(invisible(NULL))
  }
  contrast_label <- contrast_labels[[contrast_key]]

  # Long format with 5% x-offset for scrambled
  d_long <- rbind(
    data.frame(method = d$method,
               x      = d$window_size_num * 0.95,
               y      = d$cov_real_mb,
               series = "real",
               stringsAsFactors = FALSE),
    data.frame(method = d$method,
               x      = d$window_size_num * 1.05,
               y      = d$cov_scrambled_mb,
               series = "scrambled",
               stringsAsFactors = FALSE)
  )
  d_long$key <- factor(paste0(d_long$method, " (", d_long$series, ")"),
                       levels = key_levels)

  # ---- Panel A ------------------------------------------------------------
  pA <- ggplot(d_long,
               aes(x = x, y = y,
                   colour = key, linetype = key, shape = key, group = key)) +
    geom_line(linewidth = 1.0, lineend = "round") +
    geom_point(size = 3.0, stroke = 0.9, fill = "white") +
    scale_colour_manual(values = key_cols, drop = FALSE) +
    scale_linetype_manual(values = key_lty, drop = FALSE) +
    scale_shape_manual(values = key_shp, drop = FALSE) +
    x_scale +
    labs(
      title    = sprintf("A   CpG DMRs \u2014 %s (strict)", contrast_label),
      subtitle = paste0("null model: ", null_label),
      x = "Bin/window size (bp)",
      y = "DMR genome coverage (Mb)"
    ) +
    theme_radu +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1))

  # ---- Panel B: real - scrambled (delta) ----------------------------------
  pB_data <- d
  pB <- ggplot(pB_data,
               aes(x = window_size_num, y = cov_diff_mb,
                   colour = method, shape = method, group = method)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_line(linewidth = 1.0, lineend = "round") +
    geom_point(size = 3.0, stroke = 0.9, fill = "white") +
    scale_colour_manual(values = delta_cols, drop = TRUE,
                        labels = function(x) paste0(x, " (delta)"),
                        name = NULL) +
    scale_shape_manual(values = delta_shp, drop = TRUE,
                       labels = function(x) paste0(x, " (delta)"),
                       name = NULL) +
    x_scale +
    labs(
      title    = "B   CpG DMRs difference (real \u2212 scrambled) (strict)",
      subtitle = paste0("null model: ", null_label),
      x = "Bin/window size (bp)",
      y = "Delta coverage (real \u2212 scrambled, Mb)"
    ) +
    theme_radu

  # Position panel B legend depending on data sign
  if (mean(pB_data$cov_diff_mb, na.rm = TRUE) < 0) {
    pB <- pB + theme(legend.position = c(0.02, 0.02),
                     legend.justification = c(0, 0))
  } else {
    pB <- pB + theme(legend.position = c(0.98, 0.98),
                     legend.justification = c(1, 1))
  }

  caption_text <- switch(null_label,
    "Label-swap" = paste0(
      "DMRcaller-NB excluded (count constant w.r.t. window size; shown in ratio plot). ",
      "Scrambled lines offset on x-axis for visibility. ",
      "Delta \u2248 0 is the expected result for the label-swap null."),
    "Read count permutation" = paste0(
      "Real > scrambled (positive delta) confirms genuine signal. 20-seed permutation mean used as null baseline (pooled replicates)."),
    "Stratified scramble" = paste0(
      "Stratified scramble preserves per-sample read count distribution. ",
      "Scrambled lines offset on x-axis for visibility."),
    "Scrambled lines offset on x-axis for visibility."
  )

  combined <- (pA | pB) +
    plot_annotation(
      caption = caption_text,
      theme = theme(plot.caption = element_text(face = "italic",
                                                colour = "grey40", hjust = 0))
    )

  outfile_tag <- if (HAS_CONTRAST && length(contrast_levels) > 1) {
    paste0(filename_tag, "_", contrast_key)
  } else {
    filename_tag
  }
  outfile <- file.path(PLOT_DIR, paste0("radu_panels_", outfile_tag, ".pdf"))
  ggsave(outfile, combined, width = 14, height = 5.5,
         device = grDevices::cairo_pdf)
  message("Saved: ", basename(outfile))
}

# ---- Run --------------------------------------------------------------------
for (ck in contrast_levels) {
  message("\n=== Contrast: ", ck, " (", contrast_labels[[ck]], ") ===")
  make_panels("Label-swap",             "label_swap",             ck)
  make_panels("Read count permutation", "read_count_permutation", ck)
  make_panels("Stratified scramble",    "stratified_scramble",    ck)
}

message("\nDone. PDFs in: ", PLOT_DIR)
