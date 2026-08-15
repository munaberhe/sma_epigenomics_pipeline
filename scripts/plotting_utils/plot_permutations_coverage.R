.libPaths("~/R/library")
library(ggplot2)
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

# combine all readcount CSVs
csvs <- list.files("results/dmr_benchmark_mac",
                   pattern="parameter_benchmark_readcount_ws.*\\.csv",
                   full.names=TRUE)
df <- do.call(rbind, lapply(csvs, read.csv))
df <- df[df$mode=="strict",]
df$window_size_num <- as.numeric(as.character(df$window_size))

# coverage in Mb
df$cov_real_mb     <- df$n_real         * df$window_size_num / 1e6
df$cov_scr_mean_mb <- df$mean_scrambled * df$window_size_num / 1e6
df$cov_scr_sd_mb   <- df$sd_scrambled   * df$window_size_num / 1e6

cat("window sizes:", paste(sort(df$window_size_num), collapse=", "), "\n")
cat("cov_real_mb range:", range(df$cov_real_mb), "\n")

NAVY <- "#1F3A5F"; NAVY_L <- "#7FB0D3"; RED <- "#C0392B"

p <- ggplot(df, aes(x=window_size_num)) +
  geom_ribbon(aes(ymin=cov_scr_mean_mb - cov_scr_sd_mb,
                  ymax=cov_scr_mean_mb + cov_scr_sd_mb),
              fill=NAVY_L, alpha=0.35) +
  geom_line(aes(y=cov_scr_mean_mb), colour=NAVY_L,
            linewidth=1, linetype="dashed") +
  geom_point(aes(y=cov_scr_mean_mb), colour=NAVY_L, size=3) +
  geom_line(aes(y=cov_real_mb), colour=NAVY, linewidth=1.2) +
  geom_point(aes(y=cov_real_mb), colour=NAVY, size=3) +
  geom_vline(xintercept=300, linetype="dashed",
             colour=RED, linewidth=0.8) +
  annotate("text", x=320, y=max(df$cov_real_mb)*0.92,
           label="300 bp\n(selected)", colour=RED,
           hjust=0, size=3.5, fontface="bold") +
  scale_x_log10(breaks=c(100,200,300,500,1000,2000),
                labels=c("100","200","300","500","1000","2000")) +
  scale_y_continuous(labels=scales::comma) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"),
        panel.grid.minor=element_blank()) +
  labs(title="DMR genome coverage vs bin size (bins method, strict, read-count permutation)",
       subtitle=paste0("Navy = real data | dashed+band = 20-seed permutation null (mean ± SD) | ",
                       "red = selected 300 bp"),
       x="Bin size (bp)", y="DMR genome coverage (Mb)")

dir.create("results/dmr_benchmark/plots", recursive=TRUE, showWarnings=FALSE)
ggsave("results/dmr_benchmark/plots/permutations_coverage_strict.pdf",
       p, width=10, height=5.5, device=cairo_pdf)
message("saved: permutations_coverage_strict.pdf")
