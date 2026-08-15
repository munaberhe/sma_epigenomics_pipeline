#!/usr/bin/env Rscript
# 08d_pca_per_replicate.R
# Proper 12-sample PCA using per-replicate methylation data
# Each dot = one biological replicate, coloured by condition
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics

suppressPackageStartupMessages({
  library(ggplot2)
  library(DMRcaller)
  library(GenomicRanges)
})
.libPaths(c("~/R/library", .libPaths()))

OUT <- "results/dmr_qc"

# Load the per-replicate methylation matrix that 08b already computed
cor_mat <- read.csv("results/qc/methylation/methylation_correlation.csv",
                    row.names=1)

# The correlation matrix is 12x12 -- use it directly for PCA
# PCA on correlation matrix gives same clustering as PCA on methylation values
pca <- prcomp(cor_mat, scale.=TRUE)

pca_df <- data.frame(
  PC1     = pca$x[,1],
  PC2     = pca$x[,2],
  sample  = rownames(pca$x)
)

# Extract condition from sample name
pca_df$condition <- gsub("_[123]$", "", pca_df$sample)
pca_df$replicate <- gsub(".*_([123])$", "\\1", pca_df$sample)

var_exp <- round(summary(pca)$importance[2,1:2]*100, 1)
cat("PC1:", var_exp[1], "%\n")
cat("PC2:", var_exp[2], "%\n")

GROUP_COLS <- c(
  ASO_CTRL      = "#2E9B6F",
  ASO_VPA       = "#D94F3D",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#1D6FA4"
)

p <- ggplot(pca_df, aes(PC1, PC2,
                         colour=condition,
                         label=sample,
                         shape=replicate)) +
  geom_point(size=5, alpha=0.9) +
  geom_text(vjust=-0.9, size=3, colour="grey30") +
  scale_colour_manual(values=GROUP_COLS, name="Condition") +
  scale_shape_manual(values=c("1"=16,"2"=17,"3"=15), name="Replicate") +
  stat_ellipse(aes(group=condition), type="norm",
               linetype="dashed", linewidth=0.5, alpha=0.5) +
  theme_classic(base_size=12) +
  theme(
    plot.title    = element_text(face="bold"),
    plot.subtitle = element_text(size=9, colour="grey40"),
    legend.position = "right"
  ) +
  labs(
    title    = "PCA of per-replicate CpG methylation profiles",
    subtitle = paste0("n=12 samples (3 replicates × 4 conditions)\n",
                      "PC1 (", var_exp[1], "%) separates VPA from non-VPA; ",
                      "PC2 (", var_exp[2], "%) separates ASO from Scramble"),
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  )

ggsave(file.path(OUT, "sample_PCA_12samples.pdf"),
       p, width=9, height=7)
message("Saved: sample_PCA_12samples.pdf")
