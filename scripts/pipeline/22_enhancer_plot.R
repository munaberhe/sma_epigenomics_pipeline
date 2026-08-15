#!/usr/bin/env Rscript
# 22_enhancer_plot.R
# Bar plot of H9 enhancer fold enrichment across four pairwise contrasts.
# Results from 21_enhancer_enrichment_pairwise.R (50k DMR subsample, 1000 perms).

.libPaths(c("~/R/library", .libPaths()))
library(ggplot2)
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/enhancer_pairwise"

COND_COLS <- c(
  "ASO alone"  = "#1F3A5F",
  "VPA alone"  = "#F0A500",
  "ASO in VPA" = "#C0392B",
  "VPA in ASO" = "#6B4E9E"
)

df <- data.frame(
  contrast = c("ASO alone", "VPA alone", "ASO in VPA", "VPA in ASO"),
  fold     = c(1.910, 1.678, 2.332, 1.628),
  pvalue   = c(0.001, 0.001, 0.001, 0.001),
  zscore   = c(18.23, 54.15, 70.51, 49.27)
)

df$contrast <- factor(df$contrast,
  levels=c("ASO alone","VPA alone","ASO in VPA","VPA in ASO"))
df$sig <- "p < 0.001"

p <- ggplot(df, aes(x=contrast, y=fold, fill=contrast)) +
  geom_col(width=0.6) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey40", linewidth=0.5) +
  geom_text(aes(label=sprintf("%.2fx\np<0.001\n(z=%.1f)", fold, zscore)),
            vjust=-0.3, size=3.5, fontface="bold") +
  scale_fill_manual(values=COND_COLS) +
  scale_y_continuous(limits=c(0, 2.8),
                     breaks=seq(0, 2.5, 0.5)) +
  labs(x=NULL,
       y="Fold enrichment over random (H9 enhancers)",
       caption="regioneR permutation test, 1000 permutations, 50,000 DMRs per contrast (seed=42).\nH9 predicted non-promoter enhancers. Dashed line = no enrichment.") +
  theme_classic(base_size=13) +
  theme(legend.position  = "none",
        axis.text.x      = element_text(size=11),
        plot.caption     = element_text(size=8, colour="grey40"))

ggsave(file.path(OUT, "enhancer_enrichment_pairwise.pdf"),
       p, width=8, height=5, device=cairo_pdf)
message("Saved: enhancer_enrichment_pairwise.pdf")
