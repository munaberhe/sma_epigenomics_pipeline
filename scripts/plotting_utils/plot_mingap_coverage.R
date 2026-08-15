.libPaths("~/R/library")
library(ggplot2)
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

df <- read.csv("results/dmr_benchmark_mac/parameter_benchmark_neighbourhood_mingap.csv")
df <- df[df$mode=="strict" & df$window_size > 0,]
df$cov_real_mb    <- df$n_real      * 300 / 1e6
df$cov_scrambled_mb <- df$n_scrambled * 300 / 1e6

NAVY <- "#1F3A5F"; NAVY_L <- "#7FB0D3"

p <- ggplot(df, aes(x=window_size)) +
  geom_ribbon(aes(ymin=0, ymax=cov_scrambled_mb),
              fill=NAVY_L, alpha=0.4) +
  geom_line(aes(y=cov_real_mb), colour=NAVY, linewidth=1.2) +
  geom_point(aes(y=cov_real_mb), colour=NAVY, size=3) +
  geom_line(aes(y=cov_scrambled_mb), colour=NAVY_L,
            linewidth=0.8, linetype="dashed") +
  annotate("vline", xintercept=300, linetype="dashed",
           colour="#C0392B", linewidth=0.8) +
  annotate("text", x=320, y=max(df$cov_real_mb)*0.9,
           label="300 bp\n(selected)", colour="#C0392B",
           hjust=0, size=3.5, fontface="bold") +
  scale_x_continuous(breaks=c(0,100,200,300,500,1000,2000)) +
  scale_y_continuous(labels=scales::comma) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"),
        panel.grid.minor=element_blank()) +
  labs(title="DMR genome coverage vs minGap (neighbourhood method, strict)",
       subtitle="Blue line = real data | shaded band = label-swap null | red line = selected 300bp",
       x="minGap (bp)", y="DMR genome coverage (Mb)")

dir.create("results/dmr_benchmark/plots", recursive=TRUE, showWarnings=FALSE)
ggsave("results/dmr_benchmark/plots/mingap_coverage_strict.pdf",
       p, width=10, height=5.5, device=cairo_pdf)
message("saved: mingap_coverage_strict.pdf")
