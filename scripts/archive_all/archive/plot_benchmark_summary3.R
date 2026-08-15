.libPaths("~/R/library")
library(ggplot2)
library(patchwork)

# =============================================================================
# plot_benchmark_summary.R  (polished — style-matched to slides 11-15)
#
# Drop-in replacement for the previous plot_benchmark_summary.R.
# Reads the same CSVs from the same folder, writes PDFs with the same names
# into the same plots/ subfolder, and writes the same parameter-selection
# CSV table.
#
# Outputs:
#   ratio_label_swap.pdf
#   ratio_read_count_permutation.pdf
#   ratio_stratified_scramble.pdf
#   threshold_comparison_label_swap.pdf
#   threshold_comparison_read_count_permutation.pdf
#   threshold_comparison_stratified_scramble.pdf
#   benchmark_parameter_selection_table.csv  (in OUT_DIR, not PLOT_DIR)
#
# Key style fixes vs the previous version:
#   1. Bold "A" / "B" panel-title prefix where applicable
#   2. Italic grey caption footer below each plot
#   3. Navy + brick-red palette (replaces teal/blue/pink)
#   4. Filled markers with white centres (matches slides 11-15)
#   5. Log-scaled x-axis so 100/200/300 don't squash together
#   6. "Null baseline (ratio = 1)" annotation on dashed reference line
#   7. Inline boxed legend (top-right) instead of right-side legend
# =============================================================================

# ---- Where to read from ------------------------------------------------------
USE_TURBO <- as.logical(as.integer(Sys.getenv("TURBO", unset = "1")))
OUT_DIR   <- if (USE_TURBO) "results/dmr_benchmark_turbo" else "results/dmr_benchmark"
PLOT_DIR  <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Reading CSVs from: ", OUT_DIR)

csv_files <- list.files(OUT_DIR, pattern = "parameter_benchmark.*\\.csv",
                        full.names = TRUE)
if (length(csv_files) == 0) stop("No parameter_benchmark*.csv found in ", OUT_DIR)

df <- do.call(rbind, lapply(csv_files, read.csv))
df$ratio[is.infinite(df$ratio)] <- NA
df$window_size_num <- as.numeric(as.character(df$window_size))

df$scramble_method <- factor(df$scramble_method,
  levels = c("label_swap", "archie_scramble", "stratified_scramble"),
  labels = c("Label-swap", "Read count permutation", "Stratified scramble"))

df$method <- factor(df$method,
  levels = c("bins", "neighbourhood", "noise_filter"),
  labels = c("DMRcaller-B", "DMRcaller-NB", "DMRcaller-NF"))

df$mode <- factor(df$mode, levels = c("loose", "strict"))

# Gaussian kernel only for noise_filter
df <- df[!(df$method == "DMRcaller-NF" &
           df$kernel %in% c("uniform", "epanechnicov")), ]

# ---- Palette (slides 11-15) -------------------------------------------------
colours_method <- c(
  "DMRcaller-B"  = "#2C5F8D",   # navy
  "DMRcaller-NB" = "#7F8C8D",   # neutral grey
  "DMRcaller-NF" = "#C0392B")   # brick red

shapes_method <- c(
  "DMRcaller-B"  = 15,   # filled square
  "DMRcaller-NB" = 17,   # filled triangle
  "DMRcaller-NF" = 18)   # filled diamond

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

# ---- 1. Per-null signal/noise ratio plot ------------------------------------
make_ratio_plot <- function(null_label, filename_tag) {
  d <- df[df$scramble_method == null_label & df$mode == "strict", ]
  d <- d[!is.na(d$ratio), ]
  if (nrow(d) == 0) {
    message("  [skip] no rows for ", null_label)
    return(invisible(NULL))
  }

  p <- ggplot(d,
       aes(x = window_size_num, y = ratio,
           colour = method, shape = method, group = method)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    annotate("text", x = 100, y = 1, label = "  Null baseline (ratio = 1)",
             hjust = 0, vjust = -0.5, size = 3.3, colour = "grey45",
             fontface = "italic") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.8, stroke = 0.8, fill = "white") +
    scale_colour_manual(values = colours_method) +
    scale_shape_manual(values  = shapes_method) +
    x_scale +
    labs(
      title    = paste0("Signal/noise ratio \u2014 ", null_label, " (strict)"),
      subtitle = "Ratio > 1 = real data calls more DMRs than scrambled",
      x = "Bin/window size (bp)",
      y = "Signal/noise ratio (real / scrambled)"
    ) +
    theme_radu +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1))

  outfile <- file.path(PLOT_DIR, paste0("ratio_", filename_tag, ".pdf"))
  ggsave(outfile, p, width = 9, height = 5.5,
         device = grDevices::cairo_pdf)
  message("Saved: ", basename(outfile))
}

make_ratio_plot("Label-swap",             "label_swap")
make_ratio_plot("Read count permutation", "read_count_permutation")
make_ratio_plot("Stratified scramble",    "stratified_scramble")

# ---- 2. Loose vs strict comparison ------------------------------------------
make_threshold_plot <- function(null_label, filename_tag) {
  d <- df[df$scramble_method == null_label & !is.na(df$ratio), ]
  if (nrow(d) == 0) {
    message("  [skip] no rows for ", null_label)
    return(invisible(NULL))
  }

  p <- ggplot(d,
       aes(x = window_size_num, y = ratio,
           colour = method, shape = method,
           linetype = mode,
           group = interaction(method, mode))) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.8, stroke = 0.8, fill = "white") +
    scale_colour_manual(values = colours_method, name = NULL) +
    scale_shape_manual(values  = shapes_method,  name = NULL) +
    scale_linetype_manual(
      values = c("loose" = "dashed", "strict" = "solid"),
      name = NULL,
      labels = c("loose" = "Loose threshold",
                 "strict" = "Strict threshold")) +
    x_scale +
    labs(
      title    = paste0("Loose vs strict thresholds \u2014 ", null_label),
      subtitle = "Solid = strict, dashed = loose. Ratio > 1 = genuine signal.",
      x = "Bin/window size (bp)",
      y = "Signal/noise ratio (real / scrambled)",
      caption = "Strict = pval<0.01, minCpG=4, minDiff=0.2, minSize=50. Loose = pval<0.05, minCpG=1, minDiff=0.1, minSize=1."
    ) +
    theme_radu +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.box = "vertical",
          legend.spacing.y = unit(2, "pt"))

  outfile <- file.path(PLOT_DIR, paste0("threshold_comparison_", filename_tag, ".pdf"))
  ggsave(outfile, p, width = 9, height = 5.5,
         device = grDevices::cairo_pdf)
  message("Saved: ", basename(outfile))
}

make_threshold_plot("Label-swap",             "label_swap")
make_threshold_plot("Read count permutation", "read_count_permutation")
make_threshold_plot("Stratified scramble",    "stratified_scramble")

# ---- 3. Summary table at strict threshold (unchanged) -----------------------
df_strict <- df[df$mode == "strict", ]
df_strict$cov_real_mb <- df_strict$n_real * df_strict$window_size_num / 1e6

summary_table <- do.call(rbind, lapply(
  levels(df_strict$scramble_method), function(nm) {
    d <- df_strict[df_strict$scramble_method == nm, ]
    best <- do.call(rbind, lapply(levels(d$method), function(m) {
      dm <- d[d$method == m & !is.na(d$ratio), ]
      if (nrow(dm) == 0) return(NULL)
      dm[which.max(dm$ratio),
         c("method", "window_size_num", "n_real", "n_scrambled",
           "ratio", "cov_real_mb")]
    }))
    if (is.null(best)) return(NULL)
    best$null_model <- nm
    best
  }))

summary_table <- summary_table[, c("null_model", "method", "window_size_num",
                                    "n_real", "n_scrambled", "ratio",
                                    "cov_real_mb")]
colnames(summary_table) <- c("Null model", "Method", "Best window (bp)",
                              "n DMRs (real)", "n DMRs (scrambled)",
                              "Signal/noise ratio", "Coverage real (Mb)")
summary_table$`Signal/noise ratio` <- round(summary_table$`Signal/noise ratio`, 2)
summary_table$`Coverage real (Mb)` <- round(summary_table$`Coverage real (Mb)`, 1)

write.csv(summary_table,
          file.path(OUT_DIR, "benchmark_parameter_selection_table.csv"),
          row.names = FALSE)
message("\nParameter selection table:")
print(summary_table, row.names = FALSE)

message("\nAll plots saved to: ", PLOT_DIR)
