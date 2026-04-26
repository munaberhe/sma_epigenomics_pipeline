.libPaths("~/R/library")
# benchmark_compare.R
# Combines label-swap and Archie scrambling benchmark results
# Produces comparison plots for Radu meeting
# Run after both benchmark jobs finish
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(ggplot2)

OUT_DIR <- "results/dmr_benchmark"

# Load both results
label_file <- file.path(OUT_DIR, "parameter_benchmark_chr1.csv")
archie_file <- file.path(OUT_DIR, "parameter_benchmark_archie_scramble.csv")

if (!file.exists(label_file)) stop("Label-swap results not found: ", label_file)
if (!file.exists(archie_file)) stop("Archie results not found: ", archie_file)

label_df  <- read.csv(label_file)
archie_df <- read.csv(archie_file)

# Standardise columns
label_df$scramble_method  <- "Label swap"
archie_df$scramble_method <- "Archie (shuffle readsM)"

# Keep common columns
cols <- c("method", "window_size", "n_real", "n_scrambled", "ratio", "scramble_method")

# Handle label_df which may have extra columns (kernel, mode)
# Use triangular kernel, loose mode as primary comparison
if ("kernel" %in% names(label_df)) {
  label_df <- label_df[label_df$kernel == "triangular" &
                        label_df$mode == "loose", ]
}

label_df  <- label_df[, intersect(cols, names(label_df))]
archie_df <- archie_df[, intersect(cols, names(archie_df))]

combined <- rbind(label_df, archie_df)
combined$ratio_num <- as.numeric(ifelse(is.infinite(combined$ratio) |
                                        is.na(combined$ratio), NA,
                                        combined$ratio))
combined$window_size <- factor(combined$window_size,
                               levels=sort(unique(combined$window_size)))

write.csv(combined, file.path(OUT_DIR, "benchmark_combined.csv"), row.names=FALSE)

# ── Plot 1: Signal-to-noise ratio comparison ─────────────────────────────────
p1 <- ggplot(combined[!is.na(combined$ratio_num),],
             aes(x=window_size, y=ratio_num,
                 colour=scramble_method, group=scramble_method)) +
  geom_line(linewidth=1.1) +
  geom_point(size=3) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50", linewidth=0.8) +
  annotate("text", x=1, y=1.05, label="ratio = 1 (no discrimination)",
           hjust=0, size=3, colour="grey50") +
  facet_wrap(~method, scales="free_y") +
  scale_colour_manual(values=c("Label swap"="#02C39A",
                                "Archie (shuffle readsM)"="#F59E0B")) +
  labs(title="Signal-to-noise ratio — label swap vs Archie scrambling",
       subtitle="chr1 full chromosome, ASO_VPA vs ASO_CTRL, ratio > 1 = real > scrambled",
       x="Window/bin size (bp)", y="Ratio (real DMRs / scrambled DMRs)",
       colour="Scrambling method") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom")
ggsave(file.path(OUT_DIR, "benchmark_comparison_ratio.pdf"), p1, width=12, height=5)
message("Saved: benchmark_comparison_ratio.pdf")

# ── Plot 2: Real vs scrambled counts side by side ────────────────────────────
df_long <- reshape(combined,
  varying=c("n_real","n_scrambled"),
  v.names="n_dmrs", timevar="data_type",
  times=c("Real","Scrambled"), direction="long")
df_long$data_type <- factor(df_long$data_type, levels=c("Real","Scrambled"))

p2 <- ggplot(df_long[!is.na(df_long$n_dmrs),],
             aes(x=window_size, y=n_dmrs,
                 colour=data_type, linetype=scramble_method,
                 group=interaction(data_type, scramble_method))) +
  geom_line(linewidth=1.0) +
  geom_point(size=2.5) +
  facet_wrap(~method, scales="free_y") +
  scale_colour_manual(values=c(Real="#065A82", Scrambled="#DC2626")) +
  scale_linetype_manual(values=c("Label swap"="solid",
                                  "Archie (shuffle readsM)"="dashed")) +
  labs(title="DMR counts: real vs scrambled — both methods compared",
       subtitle="chr1 full chromosome, ASO_VPA vs ASO_CTRL",
       x="Window/bin size (bp)", y="Number of DMRs",
       colour="Data", linetype="Scrambling method") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom")
ggsave(file.path(OUT_DIR, "benchmark_comparison_counts.pdf"), p2, width=12, height=5)
message("Saved: benchmark_comparison_counts.pdf")

# ── Plot 3: Best parameters summary table ────────────────────────────────────
# Find best window size per method per scrambling approach
best <- aggregate(ratio_num ~ method + scramble_method,
                  data=combined[!is.na(combined$ratio_num),],
                  FUN=function(x) {
                    idx <- which.max(x)
                    x[idx]
                  })

best_ws <- do.call(rbind, lapply(split(combined[!is.na(combined$ratio_num),],
  paste(combined[!is.na(combined$ratio_num),"method"],
        combined[!is.na(combined$ratio_num),"scramble_method"])), function(df) {
  idx <- which.max(df$ratio_num)
  df[idx, c("method","scramble_method","window_size","n_real","n_scrambled","ratio_num")]
}))

message("\n=== BEST PARAMETERS BY METHOD ===")
message(sprintf("%-15s %-25s %-12s %-12s %-12s %-10s",
  "Method", "Scramble method", "Window(bp)", "Real DMRs", "Scr DMRs", "S/N Ratio"))
message(paste(rep("-", 90), collapse=""))
for (i in seq_len(nrow(best_ws))) {
  r <- best_ws[i,]
  message(sprintf("%-15s %-25s %-12s %-12s %-12s %-10s",
    r$method, r$scramble_method, r$window_size,
    r$n_real, r$n_scrambled, round(r$ratio_num, 2)))
}

write.csv(best_ws, file.path(OUT_DIR, "benchmark_best_parameters.csv"), row.names=FALSE)
message("\nSaved: benchmark_best_parameters.csv")
message("All comparison plots saved to: ", OUT_DIR)
