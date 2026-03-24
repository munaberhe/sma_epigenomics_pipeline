# ── DESeq2 Differential Expression Analysis ───────────────────────────────
# ASO1+VPA vs ASO1 — SMA Epigenomics Project

library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)

# ── 1. Setup ──────────────────────────────────────────────────────────────
setwd("~/sma_epigenomics_pipeline")

counts_file <- "results/counts/counts.txt"
padj_thresh <- 0.05
lfc_thresh  <- 1

# ── 2. Load count matrix ──────────────────────────────────────────────────
counts <- read.table(counts_file, header = TRUE, skip = 1, row.names = 1)
counts <- counts[, 6:ncol(counts)]  # drop annotation columns

head(counts)
dim(counts)

# ── 3. Sample metadata ────────────────────────────────────────────────────
# Update numbers to match your actual replicates when data arrives
coldata <- data.frame(
  condition = factor(c(rep("ASO1", ncol(counts)/2),
                       rep("ASO1_VPA", ncol(counts)/2))),
  row.names = colnames(counts)
)
coldata

# ── 4. Build DESeq2 object & run ──────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ condition)
dds <- DESeq(dds)

# ── 5. Extract results ────────────────────────────────────────────────────
res <- results(dds,
               contrast     = c("condition", "ASO1_VPA", "ASO1"),
               alpha        = padj_thresh,
               lfcThreshold = lfc_thresh)
summary(res)

res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj), ]
head(res_df)

write.csv(res_df, "results/differential/deseq2_results.csv", row.names = TRUE)

# ── 6. Volcano plot ───────────────────────────────────────────────────────
EnhancedVolcano(res_df,
                lab      = rownames(res_df),
                x        = "log2FoldChange",
                y        = "padj",
                pCutoff  = padj_thresh,
                FCcutoff = lfc_thresh,
                title    = "ASO1+VPA vs ASO1",
                subtitle = "Differential Expression")

ggsave("results/figures/volcano_plot.png", width = 10, height = 8, dpi = 150)

# ── 7. MA plot ────────────────────────────────────────────────────────────
plotMA(res,
       main  = "MA Plot — ASO1+VPA vs ASO1",
       ylim  = c(-5, 5),
       alpha = padj_thresh)

# ── 8. PCA plot ───────────────────────────────────────────────────────────
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 5) +
  geom_text(aes(label = name), vjust = -1, size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  scale_color_manual(values = c("ASO1" = "#0D1B2A", "ASO1_VPA" = "#00897B")) +
  theme_bw() +
  ggtitle("PCA — Sample Clustering by Condition")

ggsave("results/figures/pca_plot.png", width = 10, height = 8, dpi = 150)

# ── 9. Top 50 DEG heatmap ─────────────────────────────────────────────────
sig_genes <- rownames(res_df)[!is.na(res_df$padj) &
                                res_df$padj < padj_thresh][1:50]
mat <- assay(vsd)[sig_genes, ]
mat <- mat - rowMeans(mat)

annotation_col <- data.frame(
  Condition = coldata$condition,
  row.names = rownames(coldata)
)
ann_colors <- list(
  Condition = c(ASO1 = "#0D1B2A", ASO1_VPA = "#00897B")
)

pheatmap(mat,
         annotation_col    = annotation_col,
         annotation_colors = ann_colors,
         color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         fontsize_row      = 7,
         main              = "Top 50 DEGs — ASO1+VPA vs ASO1",
         filename          = "results/figures/heatmap_top50.png",
         width             = 10,
         height            = 12)

# ── 10. Sample distance heatmap ───────────────────────────────────────────
sampleDists       <- dist(t(assay(vsd)))
sampleDistMatrix  <- as.matrix(sampleDists)
colors            <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col      = colors,
         main     = "Sample-to-Sample Distances",
         filename = "results/figures/sample_distances.png",
         width    = 8,
         height   = 7)

cat("DESeq2 analysis complete. All plots saved to results/figures/\n")