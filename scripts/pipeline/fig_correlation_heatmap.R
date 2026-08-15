#!/usr/bin/env Rscript
# fig_correlation_heatmap.R
# Inter-sample CpG methylation correlation heatmap
# Scale fixed to 0.88-1.00 to show all off-diagonal values

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/thesis_figures"

# load correlation matrix
cor_mat <- read.csv("results/qc/methylation/methylation_correlation.csv",
                    row.names=1, check.names=FALSE)
cat("Correlation matrix dimensions:", dim(cor_mat), "\n")
print(round(cor_mat, 4))

# melt to long format
df <- melt(as.matrix(cor_mat), varnames=c("Sample1","Sample2"), value.name="r")

# sample order - group by condition
sample_order <- c(
  "Scramble_CTRL_1","Scramble_CTRL_2","Scramble_CTRL_3",
  "ASO_CTRL_1","ASO_CTRL_2","ASO_CTRL_3",
  "Scramble_VPA_1","Scramble_VPA_2","Scramble_VPA_3",
  "ASO_VPA_1","ASO_VPA_2","ASO_VPA_3"
)

# use whatever samples exist
existing <- intersect(sample_order, unique(df$Sample1))
if (length(existing) < 2) existing <- unique(df$Sample1)

df$Sample1 <- factor(df$Sample1, levels=existing)
df$Sample2 <- factor(df$Sample2, levels=rev(existing))

p <- ggplot(df, aes(x=Sample1, y=Sample2, fill=r)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(aes(label=sprintf("%.4f", r)), size=2.5) +
  scale_fill_gradient2(
    low="#C0392B", mid="white", high="#1F3A5F",
    midpoint=0.915,
    limits=c(0.83, 1.00),
    name="Pearson r",
    oob=scales::squish
  ) +
  labs(x=NULL, y=NULL) +
  theme_classic(base_size=11) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=8),
    axis.text.y = element_text(size=8),
    legend.position = "right"
  )

ggsave(file.path(OUT, "Fig5.1c_correlation_heatmap.pdf"),
       p, width=8, height=7, device=cairo_pdf)
message("Saved: Fig5.1c_correlation_heatmap.pdf")
