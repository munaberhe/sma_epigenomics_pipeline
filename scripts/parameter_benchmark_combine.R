.libPaths("~/R/library")
# parameter_benchmark_combine.R
# Combines bins, neighbourhood and noise_filter benchmark results
# and regenerates the final plots for the presentation
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(ggplot2)

OUT_DIR <- "results/dmr_benchmark"

bins_df  <- read.csv(file.path(OUT_DIR, "parameter_benchmark_chr1_final.csv"))
neigh_df <- read.csv(file.path(OUT_DIR, "parameter_benchmark_neighbourhood.csv"))
nf_df    <- read.csv(file.path(OUT_DIR, "parameter_benchmark_noisefilter.csv"))

summary_df <- rbind(bins_df, neigh_df, nf_df)

write.csv(summary_df,
          file.path(OUT_DIR, "parameter_benchmark_combined.csv"),
          row.names = FALSE)
message("Combined: ", nrow(summary_df), " rows")

# counts plot
summary_df$window_size <- factor(summary_df$window_size)
df_long <- reshape(summary_df,
                   varying   = c("n_real", "n_scrambled"),
                   v.names   = "n_dmrs",
                   timevar   = "type",
                   times     = c("Real", "Scrambled"),
                   direction = "long")

p1 <- ggplot(df_long, aes(x=window_size, y=n_dmrs,
                           colour=type, group=interaction(type, kernel))) +
  geom_line() + geom_point(size=2) +
  facet_grid(mode ~ method) +
  scale_colour_manual(values=c(Real="#02C39A", Scrambled="#F59E0B")) +
  labs(title = "DMR counts: real vs scrambled — chr1 first 10Mb",
       x = "Window/bin size (bp)", y = "Number of DMRs", colour = "") +
  theme_bw(base_size=11)

ggsave(file.path(OUT_DIR, "benchmark_counts_combined.pdf"), p1, width=12, height=6)
message("Saved counts plot")

# signal/noise ratio plot
summary_df$ratio_num <- as.numeric(
  ifelse(is.infinite(summary_df$ratio), NA, summary_df$ratio))

summary_df$method_kernel <- ifelse(
  summary_df$kernel == "NA" | is.na(summary_df$kernel),
  summary_df$method,
  paste0(summary_df$method, "_", summary_df$kernel))

p2 <- ggplot(summary_df, aes(x=window_size, y=ratio_num,
                              colour=method_kernel,
                              group=method_kernel)) +
  geom_line() + geom_point(size=3) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  facet_wrap(~mode) +
  labs(title = "Signal-to-noise ratio (real/scrambled) — chr1 first 10Mb",
       x = "Window/bin size (bp)", y = "Ratio (real/scrambled)",
       colour = "Method/Kernel") +
  theme_bw(base_size=11)

ggsave(file.path(OUT_DIR, "benchmark_ratio_combined.pdf"), p2, width=12, height=5)
message("Saved ratio plot")

message("\nDone. Combined outputs in: ", OUT_DIR)
