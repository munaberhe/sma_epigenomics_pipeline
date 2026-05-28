.libPaths("~/R/library")
library(ggplot2)
library(patchwork)
library(scales)

# =============================================================================
# plot_mingap_results.R  (polished v3 — slide 14 style)
# SMA Epigenomics Pipeline — Muna Berhe, QMUL
#
# Drop-in replacement for the previous plot_mingap_results.R.
#
# v3 fixes vs v2:
#   - Stratified scramble dropped (per user request)
#   - SEPARATE plot per null model (was one facetted plot which got squashed
#     because Read-count scrambled goes to 50,000 while Label-swap scrambled
#     stays under 1,000)
#   - Each null model plot has its own y-axis so series are readable
# =============================================================================

OUT_DIR  <- "results/dmr_benchmark"
PLOT_DIR <- "results/dmr_benchmark/plots_v2"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap_final.csv"))

df$mingap_num <- as.numeric(as.character(df$mingap))
df$minsize    <- factor(df$minsize, labels = c("minSize=100", "minSize=200"))
df$scramble_method <- factor(df$scramble_method,
  levels = c("label_swap", "archie", "stratified"),
  labels = c("Label-swap", "Read count permutation", "Stratified scramble"))

# Drop Stratified
df <- df[df$scramble_method != "Stratified scramble", ]
df$scramble_method <- droplevels(df$scramble_method)

# ---- Palette per null -------------------------------------------------------
palettes <- list(
  "Label-swap" = list(
    real_col  = "#2C5F8D", scram_col = "#7FB0D3"),
  "Read count permutation" = list(
    real_col  = "#C0392B", scram_col = "#E8A39B")
)

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
    legend.background  = element_rect(fill = scales::alpha("white", 0.92),
                                      colour = "grey85", linewidth = 0.3),
    legend.margin      = margin(2, 6, 2, 6),
    plot.margin        = margin(10, 14, 10, 14)
  )

x_scale <- scale_x_log10(
  breaks = sort(unique(df$mingap_num)),
  labels = sort(unique(df$mingap_num)),
  expand = expansion(mult = c(0.05, 0.05))
)

# ---- One plot per null model ------------------------------------------------
make_plot <- function(null_name, file_suffix) {
  d <- df[df$scramble_method == null_name, ]
  if (nrow(d) == 0) return(invisible(NULL))

  p <- palettes[[null_name]]

  # Long format + 5% offset
  d_long <- rbind(
    data.frame(mingap_num = d$mingap_num, minsize = d$minsize,
               x = d$mingap_num * 0.95, y = d$n_real,
               series = paste0(null_name, " (real)"),
               stringsAsFactors = FALSE),
    data.frame(mingap_num = d$mingap_num, minsize = d$minsize,
               x = d$mingap_num * 1.05, y = d$n_scrambled,
               series = paste0(null_name, " (scrambled)"),
               stringsAsFactors = FALSE)
  )
  d_long$series <- factor(d_long$series,
    levels = c(paste0(null_name, " (real)"),
               paste0(null_name, " (scrambled)")))

  series_cols <- setNames(c(p$real_col, p$scram_col), levels(d_long$series))
  series_lty  <- setNames(c("solid", "longdash"),     levels(d_long$series))
  series_shp  <- setNames(c(15, 0),                    levels(d_long$series))

  pl <- ggplot(d_long,
         aes(x = x, y = y,
             colour = series, linetype = series, shape = series,
             group = interaction(series, minsize))) +
    geom_line(linewidth = 1.0, lineend = "round") +
    geom_point(size = 2.8, stroke = 0.9, fill = "white") +
    facet_wrap(~minsize, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = series_cols) +
    scale_linetype_manual(values = series_lty) +
    scale_shape_manual(values = series_shp) +
    x_scale +
    labs(
      title    = sprintf("Neighbourhood method \u2014 Real vs Scrambled DMR Counts (%s)",
                          null_name),
      subtitle = "strict threshold, minDiff=0.2",
      x = "minGap (bp)",
      y = "Number of DMRs",
      caption = "DMRcaller neighbourhood method. Strict threshold, minDiff=0.2. Scrambled lines offset on x-axis for visibility."
    ) +
    theme_radu +
    theme(legend.position = "top")

  outfile <- file.path(PLOT_DIR,
    sprintf("mingap_benchmark_counts_%s.pdf", file_suffix))
  ggsave(outfile, pl, width = 12, height = 6,
         device = grDevices::cairo_pdf)
  message("Saved: ", basename(outfile))
}

make_plot("Label-swap",             "label_swap")
make_plot("Read count permutation", "read_count_permutation")
