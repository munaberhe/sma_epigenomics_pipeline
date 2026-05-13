.libPaths("~/R/library")
# plot_mingap_results.R
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(ggplot2)
library(patchwork)

OUT_DIR  <- "results/dmr_benchmark"
PLOT_DIR <- "results/dmr_benchmark/plots_v2"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap_final.csv"))
df$mingap  <- factor(df$mingap,  levels = sort(unique(df$mingap)))
df$minsize <- factor(df$minsize, labels = c("minSize=100", "minSize=200"))
df$scramble_method <- factor(df$scramble_method,
  levels = c("label_swap", "archie", "stratified"),
  labels = c("Label-swap", "Read count permutation", "Stratified scramble"))

colours_null <- c(
  "Label-swap"             = "#1565C0",
  "Read count permutation" = "#AD1457",
  "Stratified scramble"    = "#2E7D32")

theme_pub <- theme_bw(base_size = 13) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text       = element_text(face = "bold"))

pA <- ggplot(df,
    aes(x = mingap, y = n_real, colour = scramble_method,
        shape = scramble_method, group = scramble_method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~minsize, ncol = 2) +
  scale_colour_manual(values = colours_null, name = "Null model") +
  scale_shape_manual(values = c(15, 17, 18), name = "Null model") +
  labs(title = "A   Neighbourhood method — real DMR counts",
       x = "minGap (bp)", y = "n DMRs (real)") +
  theme_pub

pB <- ggplot(df,
    aes(x = mingap, y = n_scrambled, colour = scramble_method,
        shape = scramble_method, group = scramble_method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~minsize, ncol = 2) +
  scale_colour_manual(values = colours_null, name = "Null model") +
  scale_shape_manual(values = c(15, 17, 18), name = "Null model") +
  labs(title = "B   Neighbourhood method — scrambled DMR counts",
       x = "minGap (bp)", y = "n DMRs (scrambled)") +
  theme_pub

combined <- pA / pB + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(PLOT_DIR, "mingap_benchmark_counts.pdf"),
       combined, width = 10, height = 8)
message("Saved: mingap_benchmark_counts.pdf")
