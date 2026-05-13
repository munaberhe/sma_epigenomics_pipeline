.libPaths("~/R/library")
library(ggplot2)
library(patchwork)

OUT_DIR  <- "results/dmr_benchmark"
PLOT_DIR <- "results/dmr_benchmark/plots"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

csv_files <- list.files(OUT_DIR, pattern = "parameter_benchmark.*\\.csv", full.names = TRUE)
df <- do.call(rbind, lapply(csv_files, read.csv))
df$ratio[is.infinite(df$ratio)] <- NA
df$window_size_num <- as.numeric(as.character(df$window_size))

df$scramble_method <- factor(df$scramble_method,
  levels = c("label_swap","archie_scramble","stratified_scramble"),
  labels = c("Label-swap","Read count permutation","Stratified scramble"))

df$method <- factor(df$method,
  levels = c("bins","neighbourhood","noise_filter"),
  labels = c("DMRcaller-B","DMRcaller-NB","DMRcaller-NF"))

df$mode <- factor(df$mode, levels = c("loose","strict"))

# gaussian kernel only for noise_filter
df <- df[!(df$method == "DMRcaller-NF" &
           df$kernel %in% c("uniform","epanechnicov")), ]

colours_method <- c(
  "DMRcaller-B"  = "#00897B",
  "DMRcaller-NB" = "#1565C0",
  "DMRcaller-NF" = "#AD1457")

shapes_method <- c(
  "DMRcaller-B"  = 0,
  "DMRcaller-NB" = 2,
  "DMRcaller-NF" = 5)

theme_pub <- theme_bw(base_size = 13) +
  theme(legend.position  = "right",
        panel.grid.minor = element_blank())

x_scale <- scale_x_continuous(
  breaks = c(100, 200, 300, 500, 1000, 2000),
  labels = c("100","200","300","500","1000","2000"))

# ── 1. Signal/noise ratio — one PDF per null model ────────────────────────────
make_ratio_plot <- function(null_label, filename_tag) {
  d <- df[df$scramble_method == null_label & df$mode == "strict", ]
  d <- d[!is.na(d$ratio), ]

  p <- ggplot(d,
      aes(x = window_size_num, y = ratio,
          colour = method, shape = method, group = method)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    x_scale +
    scale_colour_manual(values = colours_method, name = "Method") +
    scale_shape_manual(values  = shapes_method,  name = "Method") +
    labs(title    = paste0("Signal/noise ratio — ", null_label, " (strict)"),
         subtitle = "Ratio > 1 = real data calls more DMRs than scrambled",
         x = "bin/window size (bp)",
         y = "Signal/noise ratio (real / scrambled)") +
    theme_pub

  outfile <- file.path(PLOT_DIR, paste0("ratio_", filename_tag, ".pdf"))
  ggsave(outfile, p, width = 8, height = 5)
  message("Saved: ", basename(outfile))
}

make_ratio_plot("Label-swap",             "label_swap")
make_ratio_plot("Read count permutation", "read_count_permutation")
make_ratio_plot("Stratified scramble",    "stratified_scramble")

# ── 2. Loose vs strict comparison — one plot per null model ───────────────────
make_threshold_plot <- function(null_label, filename_tag) {
  d <- df[df$scramble_method == null_label & !is.na(df$ratio), ]

  p <- ggplot(d,
      aes(x = window_size_num, y = ratio,
          colour = method, shape = method,
          linetype = mode,
          group = interaction(method, mode))) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    x_scale +
    scale_colour_manual(values = colours_method, name = "Method") +
    scale_shape_manual(values  = shapes_method,  name = "Method") +
    scale_linetype_manual(values = c("loose" = "dashed", "strict" = "solid"),
                          name = "Threshold") +
    labs(title    = paste0("Loose vs strict thresholds — ", null_label),
         subtitle = "Solid = strict, dashed = loose. Ratio > 1 = genuine signal.",
         x = "bin/window size (bp)",
         y = "Signal/noise ratio (real / scrambled)") +
    theme_pub

  outfile <- file.path(PLOT_DIR, paste0("threshold_comparison_", filename_tag, ".pdf"))
  ggsave(outfile, p, width = 8, height = 5)
  message("Saved: ", basename(outfile))
}

make_threshold_plot("Label-swap",             "label_swap")
make_threshold_plot("Read count permutation", "read_count_permutation")
make_threshold_plot("Stratified scramble",    "stratified_scramble")

# ── 3. Summary table at strict threshold ─────────────────────────────────────
df_strict <- df[df$mode == "strict", ]
df_strict$cov_real_mb <- df_strict$n_real * df_strict$window_size_num / 1e6

summary_table <- do.call(rbind, lapply(
  levels(df_strict$scramble_method), function(nm) {
    d <- df_strict[df_strict$scramble_method == nm, ]
    # find window size with max ratio per method
    best <- do.call(rbind, lapply(levels(d$method), function(m) {
      dm <- d[d$method == m & !is.na(d$ratio), ]
      if (nrow(dm) == 0) return(NULL)
      dm[which.max(dm$ratio), c("method","window_size_num","n_real","n_scrambled","ratio","cov_real_mb")]
    }))
    if (is.null(best)) return(NULL)
    best$null_model <- nm
    best
  }))

summary_table <- summary_table[, c("null_model","method","window_size_num",
                                    "n_real","n_scrambled","ratio","cov_real_mb")]
colnames(summary_table) <- c("Null model","Method","Best window (bp)",
                              "n DMRs (real)","n DMRs (scrambled)",
                              "Signal/noise ratio","Coverage real (Mb)")
summary_table$`Signal/noise ratio` <- round(summary_table$`Signal/noise ratio`, 2)
summary_table$`Coverage real (Mb)` <- round(summary_table$`Coverage real (Mb)`, 1)

write.csv(summary_table,
          file.path(OUT_DIR, "benchmark_parameter_selection_table.csv"),
          row.names = FALSE)
message("\nParameter selection table:")
print(summary_table, row.names = FALSE)

message("\nAll plots saved to: ", PLOT_DIR)
