# ── DMR-DEG Integration ───────────────────────────────────────────────────
# ASO1+VPA vs ASO1 — SMA Epigenomics Project

library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)

# ── 1. Setup ──────────────────────────────────────────────────────────────
setwd("~/sma_epigenomics_pipeline")

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ── 2. Load DMRs into GRanges ─────────────────────────────────────────────
dmrs_df <- read.csv("results/differential/dmrs.csv")

dmrs_gr <- makeGRangesFromDataFrame(dmrs_df,
                                    keep.extra.columns = TRUE,
                                    seqnames.field     = "seqnames",
                                    start.field        = "start",
                                    end.field          = "end")
dmrs_gr

# ── 3. Load significant DEGs ──────────────────────────────────────────────
degs_df  <- read.csv("results/differential/deseq2_results.csv", row.names = 1)
degs_sig <- degs_df[!is.na(degs_df$padj) & degs_df$padj < 0.05, ]
cat("Significant DEGs:", nrow(degs_sig), "\n")

# ── 4. Load gene annotation & match to DEGs ───────────────────────────────
gtf       <- import("data/reference/hg38.ensGene.gtf")
genes     <- gtf[gtf$type == "gene"]
deg_genes <- genes[genes$gene_name %in% rownames(degs_sig)]
cat("DEGs with coordinates:", length(deg_genes), "\n")

# ── 5. Find DMR-DEG overlaps ──────────────────────────────────────────────
hits <- findOverlaps(dmrs_gr, deg_genes)
cat("DMR-DEG overlaps found:", length(hits), "\n")

overlap_df <- data.frame(
  dmr_seqnames     = as.character(seqnames(dmrs_gr)[queryHits(hits)]),
  dmr_start        = start(dmrs_gr)[queryHits(hits)],
  dmr_end          = end(dmrs_gr)[queryHits(hits)],
  methylation_diff = dmrs_gr$methylationDifference[queryHits(hits)],
  gene_name        = deg_genes$gene_name[subjectHits(hits)],
  log2FoldChange   = degs_sig[deg_genes$gene_name[subjectHits(hits)], "log2FoldChange"],
  padj             = degs_sig[deg_genes$gene_name[subjectHits(hits)], "padj"]
)

head(overlap_df)
write.csv(overlap_df, "results/differential/dmr_deg_overlap.csv", row.names = FALSE)

# ── 6. DMR-DEG scatter plot ───────────────────────────────────────────────
ggplot(overlap_df, aes(x = methylation_diff, y = log2FoldChange)) +
  geom_point(aes(color = -log10(padj)), size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_gradient(low = "#B2DFDB", high = "#0D1B2A",
                       name = "-log10(padj)") +
  geom_text(data = overlap_df[abs(overlap_df$log2FoldChange) > 2, ],
            aes(label = gene_name), size = 2.5, vjust = -0.5) +
  theme_bw() +
  labs(title    = "DMR Methylation Change vs Gene Expression Change",
       subtitle = "Each point = one DMR-DEG overlap",
       x        = "Methylation Difference (ASO1+VPA - ASO1)",
       y        = "log2 Fold Change (ASO1+VPA / ASO1)")

ggsave("results/figures/dmr_deg_scatter.png", width = 10, height = 8, dpi = 150)

# ── 7. Genomic annotation of DMRs ─────────────────────────────────────────
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno <- annotatePeak(dmrs_gr,
                     tssRegion = c(-3000, 3000),
                     TxDb      = txdb,
                     annoDb    = "org.Hs.eg.db")

# Pie chart — where do DMRs fall in the genome?
plotAnnoPie(anno, main = "Genomic Distribution of DMRs")

# Bar chart version
plotAnnoBar(anno, title = "Genomic Annotation of DMRs")

# Distance to nearest TSS
plotDistToTSS(anno, title = "Distance of DMRs to Nearest TSS")

cat("Integration complete. All plots saved to results/figures/\n")