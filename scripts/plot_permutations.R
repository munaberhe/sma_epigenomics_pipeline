.libPaths("~/R/library")
library(ggplot2)
library(patchwork)

df <- read.csv("results/dmr_benchmark/benchmark_permutations_summary.csv")
PLOT_DIR <- "results/dmr_benchmark/plots"
dir.create(PLOT_DIR, showWarnings=FALSE, recursive=TRUE)

df$cov_real_mb     <- df$n_real         * df$window_size / 1e6
df$cov_scr_mean_mb <- df$mean_scrambled * df$window_size / 1e6
df$cov_scr_sd_mb   <- df$sd_scrambled   * df$window_size / 1e6
df$cov_diff_mb     <- df$cov_real_mb - df$cov_scr_mean_mb

theme_pub <- theme_bw(base_size=13) +
  theme(legend.position="right", panel.grid.minor=element_blank())

col <- "#00897B"

pA <- ggplot(df, aes(x=window_size)) +
  geom_ribbon(aes(ymin=cov_scr_mean_mb - cov_scr_sd_mb,
                  ymax=cov_scr_mean_mb + cov_scr_sd_mb),
              fill=col, alpha=0.15) +
  geom_line(aes(y=cov_scr_mean_mb), colour=col, linewidth=1, linetype="dashed") +
  geom_point(aes(y=cov_scr_mean_mb), colour=col, shape=0, size=2.5) +
  geom_line(aes(y=cov_real_mb), colour=col, linewidth=1, linetype="solid") +
  geom_point(aes(y=cov_real_mb), colour=col, shape=0, size=2.5) +
  scale_x_continuous(breaks=c(100,200,300,500,1000,2000)) +
  annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5,
           label="solid = Real\ndashed = Scrambled mean\u00b1SD", size=3.2, colour="grey40") +
  labs(title="A   CpG DMRs - ASO_VPA vs ASO_CTRL (strict)\nnull model: Read count permutation (20 seeds, DMRcaller-B)",
       x="bin/window size (bp)", y="DMR genome coverage (Mb)") +
  theme_pub

pB <- ggplot(df, aes(x=window_size)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_ribbon(aes(ymin=cov_diff_mb - cov_scr_sd_mb,
                  ymax=cov_diff_mb + cov_scr_sd_mb),
              fill=col, alpha=0.15) +
  geom_line(aes(y=cov_diff_mb), colour=col, linewidth=1) +
  geom_point(aes(y=cov_diff_mb), colour=col, shape=0, size=2.5) +
  scale_x_continuous(breaks=c(100,200,300,500,1000,2000)) +
  labs(title="B   CpG DMRs difference (real minus scrambled mean)\nnull model: Read count permutation (20 seeds, DMRcaller-B)",
       x="bin/window size (bp)", y="Delta coverage (Mb)") +
  theme_pub

combined <- pA + plot_spacer() + pB + plot_layout(widths=c(1,0.05,1))
ggsave(file.path(PLOT_DIR, "radu_panels_read_count_permutation_20seeds.pdf"),
       combined, width=14, height=5)
message("Saved: radu_panels_read_count_permutation_20seeds.pdf")
