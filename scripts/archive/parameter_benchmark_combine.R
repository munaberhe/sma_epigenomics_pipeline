.libPaths("~/R/library")
library(ggplot2)

OUT_DIR  <- "results/dmr_benchmark"
PLOT_DIR <- "results/dmr_benchmark/plots"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

csv_files <- list.files(OUT_DIR, pattern = "parameter_benchmark.*\\.csv", full.names = TRUE)
message("Found CSVs: ", paste(basename(csv_files), collapse=", "))

df <- do.call(rbind, lapply(csv_files, read.csv))
message("Total rows: ", nrow(df))

df$ratio[is.infinite(df$ratio)] <- NA
df$window_size_num <- as.numeric(as.character(df$window_size))
df$window_size     <- factor(df$window_size_num, levels=sort(unique(df$window_size_num)))

df$scramble_method <- factor(df$scramble_method,
  levels = c("label_swap","archie_scramble","stratified_scramble"),
  labels = c("Label-swap","Read count permutation","Stratified scramble"))

df$method <- factor(df$method,
  levels = c("bins","neighbourhood","noise_filter"),
  labels = c("Bins","Neighbourhood","Noise_filter"))

df$cov_real      <- df$n_real      * df$window_size_num
df$cov_scrambled <- df$n_scrambled * df$window_size_num
df$cov_diff      <- df$cov_real    - df$cov_scrambled

colours_null <- c(
  "Label-swap"             = "#1565C0",
  "Read count permutation" = "#AD1457",
  "Stratified scramble"    = "#2E7D32")

colours_method <- c(
  "Bins"         = "#00897B",
  "Neighbourhood"= "#FB8C00",
  "Noise_filter" = "#6A1B9A")

theme_pub <- theme_bw(base_size = 13) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="#0D2137"),
        strip.text       = element_text(colour="white", face="bold"))

# Strict only for clean plots
df_strict <- df[df$mode == "strict",]
# Bins only (one kernel per noise_filter)
df_strict_nf <- df_strict[!(df_strict$method=="Noise_filter" & df_strict$kernel %in% c("uniform","epanechnicov")),]

# --- RADU PANEL A: genome coverage real vs scrambled ---
# One plot per null model, methods as coloured lines
# Solid = real, dashed = scrambled
p_panelA <- ggplot(df_strict_nf,
    aes(x=window_size_num, colour=method, group=method)) +
  geom_line(aes(y=cov_real),      linewidth=1.2, linetype="solid") +
  geom_line(aes(y=cov_scrambled), linewidth=1.2, linetype="dashed") +
  geom_point(aes(y=cov_real),      size=2.5) +
  geom_point(aes(y=cov_scrambled), size=2.5, shape=1) +
  facet_wrap(~scramble_method, ncol=3) +
  scale_x_continuous(breaks=c(100,200,300,500,1000,2000),
                     trans="log10",
                     labels=c("100","200","300","500","1000","2000")) +
  scale_colour_manual(values=colours_method) +
  labs(title="Panel A — Genome coverage of DMRs (strict thresholds)",
       subtitle="Solid = real data, dashed = scrambled. Coverage approximated as n_DMRs x window_size.",
       x="Window/bin size (bp, log scale)",
       y="Approximate genome coverage (bp)",
       colour="Method") +
  theme_pub

ggsave(file.path(PLOT_DIR, "radu_panelA_coverage.pdf"),
       p_panelA, width=14, height=5)
message("Saved: radu_panelA_coverage.pdf")

# --- RADU PANEL B: coverage difference (real minus scrambled) ---
p_panelB <- ggplot(df_strict_nf,
    aes(x=window_size_num, y=cov_diff, colour=method, group=method)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_line(linewidth=1.2) +
  geom_point(size=2.5) +
  facet_wrap(~scramble_method, ncol=3) +
  scale_x_continuous(breaks=c(100,200,300,500,1000,2000),
                     trans="log10",
                     labels=c("100","200","300","500","1000","2000")) +
  scale_colour_manual(values=colours_method) +
  labs(title="Panel B — Coverage difference: real minus scrambled (strict thresholds)",
       subtitle="Above zero = real data has more DMR coverage than scrambled = genuine signal",
       x="Window/bin size (bp, log scale)",
       y="Coverage difference (bp)",
       colour="Method") +
  theme_pub

ggsave(file.path(PLOT_DIR, "radu_panelB_difference.pdf"),
       p_panelB, width=14, height=5)
message("Saved: radu_panelB_difference.pdf")

# --- SIGNAL/NOISE RATIO: all null models, all methods ---
df_ratio <- df_strict_nf[!is.na(df_strict_nf$ratio) & is.finite(df_strict_nf$ratio),]

p_ratio <- ggplot(df_ratio,
    aes(x=window_size_num, y=ratio, colour=scramble_method, group=scramble_method)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  geom_line(linewidth=1.2) +
  geom_point(size=2.5) +
  facet_wrap(~method, ncol=3, scales="free_y") +
  scale_x_continuous(breaks=c(100,200,300,500,1000,2000),
                     trans="log10",
                     labels=c("100","200","300","500","1000","2000")) +
  scale_colour_manual(values=colours_null) +
  labs(title="Signal/noise ratio — strict thresholds",
       subtitle="Ratio > 1 = real data calls more DMRs than scrambled",
       x="Window/bin size (bp, log scale)",
       y="Signal/noise ratio (real / scrambled)",
       colour="Null model") +
  theme_pub

ggsave(file.path(PLOT_DIR, "signal_noise_ratio.pdf"),
       p_ratio, width=14, height=5)
message("Saved: signal_noise_ratio.pdf")

# --- SUPPLEMENTARY: loose thresholds faceted ---
df_loose <- df[df$mode=="loose",]
df_loose_nf <- df_loose[!(df_loose$method=="Noise_filter" & df_loose$kernel %in% c("uniform","epanechnicov")),]
df_loose_ratio <- df_loose_nf[!is.na(df_loose_nf$ratio) & is.finite(df_loose_nf$ratio),]

p_supp <- ggplot(df_loose_ratio,
    aes(x=window_size_num, y=ratio, colour=scramble_method, group=scramble_method)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  geom_line(linewidth=1.2) +
  geom_point(size=2.5) +
  facet_wrap(~method, ncol=3, scales="free_y") +
  scale_x_continuous(breaks=c(100,200,300,500,1000,2000),
                     trans="log10",
                     labels=c("100","200","300","500","1000","2000")) +
  scale_colour_manual(values=colours_null) +
  labs(title="Signal/noise ratio — loose thresholds (supplementary)",
       subtitle="Ratio > 1 = real data calls more DMRs than scrambled",
       x="Window/bin size (bp, log scale)",
       y="Signal/noise ratio (real / scrambled)",
       colour="Null model") +
  theme_pub

ggsave(file.path(PLOT_DIR, "signal_noise_ratio_loose_supplementary.pdf"),
       p_supp, width=14, height=5)
message("Saved: signal_noise_ratio_loose_supplementary.pdf")

# --- Summary table ---
message("\n=== SIGNAL/NOISE SUMMARY (strict) ===")
summary_wide <- reshape(
  df_strict_nf[,c("method","window_size","scramble_method","ratio")],
  idvar="method", timevar=c("scramble_method"),
  direction="wide")
write.csv(df_strict_nf,
          file.path(OUT_DIR, "benchmark_summary_wide.csv"),
          row.names=FALSE)
message("\nAll plots saved to: ", PLOT_DIR)
