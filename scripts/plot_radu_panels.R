.libPaths("~/R/library")
library(ggplot2)
library(patchwork)

USE_TURBO <- as.logical(as.integer(Sys.getenv("TURBO", unset = "1")))
OUT_DIR   <- if (USE_TURBO) "results/dmr_benchmark_turbo" else "results/dmr_benchmark"
PLOT_DIR  <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

csv_files <- list.files(OUT_DIR, pattern = "parameter_benchmark.*\\.csv", full.names = TRUE)
df <- do.call(rbind, lapply(csv_files, read.csv))
df$ratio[is.infinite(df$ratio)] <- NA
df$window_size_num <- as.numeric(as.character(df$window_size))

df$method <- factor(df$method,
  levels = c("bins","neighbourhood","noise_filter"),
  labels = c("DMRcaller-B","DMRcaller-NB","DMRcaller-NF"))

df$scramble_method <- factor(df$scramble_method,
  levels = c("label_swap","archie_scramble","stratified_scramble"),
  labels = c("Label-swap","Read count permutation","Stratified scramble"))

# strict + gaussian kernel for noise_filter
df_s <- df[df$mode == "strict", ]
df_s <- df_s[!(df_s$method == "DMRcaller-NF" &
               df_s$kernel %in% c("uniform","epanechnicov")), ]

df_s$cov_real_mb      <- df_s$n_real      * df_s$window_size_num / 1e6
df_s$cov_scrambled_mb <- df_s$n_scrambled * df_s$window_size_num / 1e6
df_s$cov_diff_mb      <- df_s$cov_real_mb - df_s$cov_scrambled_mb

colours_method <- c(
  "DMRcaller-B"  = "#00897B",
  "DMRcaller-NB" = "#1565C0",
  "DMRcaller-NF" = "#AD1457")

shapes_method <- c(
  "DMRcaller-B"  = 0,
  "DMRcaller-NB" = 2,
  "DMRcaller-NF" = 5)

all_labels <- c("DMRcaller-B", "DMRcaller-NB", "DMRcaller-NF",
                "DMRcaller-B (scrambled)", "DMRcaller-NB (scrambled)", "DMRcaller-NF (scrambled)")

colours_all <- setNames(c(colours_method, colours_method), all_labels)
shapes_all  <- setNames(c(shapes_method,  shapes_method),  all_labels)
lines_all   <- setNames(c("solid","solid","solid","dashed","dashed","dashed"), all_labels)

theme_pub <- theme_bw(base_size = 13) +
  theme(legend.position  = "right",
        panel.grid.minor = element_blank())

x_scale <- scale_x_continuous(
  breaks = c(100, 200, 300, 500, 1000, 2000),
  labels = c("100","200","300","500","1000","2000"))

make_panels <- function(null_model_label, filename_tag) {
  d <- df_s[df_s$scramble_method == null_model_label, ]

  # Panel A — reshape long
  df_long <- rbind(
    data.frame(d[, c("method","window_size_num")],
               cov_mb = d$cov_real_mb, data_type = "Real"),
    data.frame(d[, c("method","window_size_num")],
               cov_mb = d$cov_scrambled_mb, data_type = "Scrambled")
  )
  df_long$label <- factor(
    paste0(df_long$method, ifelse(df_long$data_type == "Scrambled", " (scrambled)", "")),
    levels = all_labels)

  pA <- ggplot(df_long,
      aes(x = window_size_num, y = cov_mb,
          colour = label, shape = label, linetype = label, group = label)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    x_scale +
    scale_colour_manual(values = colours_all, name = NULL) +
    scale_shape_manual(values  = shapes_all,  name = NULL) +
    scale_linetype_manual(values = lines_all, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, " Mb")) +
    labs(title    = paste0("A   CpG DMRs — ASO_VPA vs ASO_CTRL\nnull model: ", null_model_label),
         x = "bin/window size (bp)",
         y = "DMR genome coverage (Mb)") +
    theme_pub

  pB <- ggplot(d,
      aes(x = window_size_num, y = cov_diff_mb,
          colour = method, shape = method, group = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    x_scale +
    scale_colour_manual(
      values = colours_method,
      labels = paste0(names(colours_method), " (difference)"),
      name   = NULL) +
    scale_shape_manual(
      values = shapes_method,
      labels = paste0(names(shapes_method), " (difference)"),
      name   = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, " Mb")) +
    labs(title = paste0("B   CpG DMRs difference (real minus scrambled)\nnull model: ", null_model_label),
         x = "bin/window size (bp)",
         y = "DMR genome coverage (Mb)") +
    theme_pub

  combined <- pA + pB + plot_layout(ncol = 2)

  outfile <- file.path(PLOT_DIR, paste0("radu_panels_", filename_tag, ".pdf"))
  ggsave(outfile, combined, width = 16, height = 6)
  message("Saved: ", basename(outfile))
}

make_panels("Label-swap",             "label_swap")
make_panels("Read count permutation", "read_count_permutation")
make_panels("Stratified scramble",    "stratified_scramble")

message("Done.")
