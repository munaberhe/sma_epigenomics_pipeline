# DMR-DEG integration and pathway enrichment
# SMA epigenomics project - Muna Berhe
# Overlapping differential methylation and expression results

library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(tidyverse)

setwd("~/sma_epigenomics_pipeline")

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/enrichment", showWarnings = FALSE, recursive = TRUE)

# load DMRs and convert to GRanges object
dmrs_df <- read.csv("results/differential/dmrs.csv")

dmrs_gr <- makeGRangesFromDataFrame(dmrs_df,
                                    keep.extra.columns = TRUE,
                                    seqnames.field     = "seqnames",
                                    start.field        = "start",
                                    end.field          = "end")
cat("DMRs loaded:", length(dmrs_gr), "\n")

# load DESeq2 results and filter for significance
degs_df  <- read.csv("results/differential/deseq2_results.csv", row.names = 1)
degs_sig <- degs_df[!is.na(degs_df$padj) & degs_df$padj < 0.05, ]
cat("Significant DEGs:", nrow(degs_sig), "\n")

# get genomic coordinates for significant DEGs from the GTF
gtf       <- import("data/reference/hg38.ensGene.gtf")
genes     <- gtf[gtf$type == "gene"]
deg_genes <- genes[genes$gene_name %in% rownames(degs_sig)]
cat("DEGs with coordinates:", length(deg_genes), "\n")

# find overlaps between DMRs and DEG loci
hits <- findOverlaps(dmrs_gr, deg_genes)
cat("DMR-DEG overlaps:", length(hits), "\n")

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

# scatter plot of methylation change vs expression change
# negative correlation expected if methylation is driving expression changes
ggplot(overlap_df, aes(x = methylation_diff, y = log2FoldChange)) +
  geom_point(aes(color = -log10(padj)), size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_gradient(low = "#B2DFDB", high = "#0D1B2A",
                       name = "-log10(padj)") +
  geom_text(data = overlap_df[abs(overlap_df$log2FoldChange) > 2, ],
            aes(label = gene_name), size = 2.5, vjust = -0.5) +
  theme_bw() +
  labs(title    = "Methylation change vs expression change",
       subtitle = "Each point = one DMR-DEG overlap",
       x        = "Methylation difference (ASO1+VPA - ASO1)",
       y        = "log2 fold change (ASO1+VPA / ASO1)")

ggsave("results/figures/dmr_deg_scatter.png", width = 10, height = 8, dpi = 150)

# annotate DMRs to genomic features using ChIPseeker
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno <- annotatePeak(dmrs_gr,
                     tssRegion = c(-3000, 3000),
                     TxDb      = txdb,
                     annoDb    = "org.Hs.eg.db")

plotAnnoPie(anno, main = "Genomic distribution of DMRs")
ggsave("results/figures/dmr_annotation_pie.png", width = 8, height = 8, dpi = 150)

plotAnnoBar(anno, title = "Genomic annotation of DMRs")
ggsave("results/figures/dmr_annotation_bar.png", width = 10, height = 6, dpi = 150)

plotDistToTSS(anno, title = "Distance of DMRs to nearest TSS")
ggsave("results/figures/dmr_distance_tss.png", width = 10, height = 6, dpi = 150)

# convert gene symbols to Entrez IDs for clusterProfiler
deg_entrez <- bitr(rownames(degs_sig),
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)
cat("DEGs mapped to Entrez IDs:", nrow(deg_entrez), "\n")

dmr_gene_entrez <- bitr(unique(overlap_df$gene_name),
                         fromType = "SYMBOL",
                         toType   = "ENTREZID",
                         OrgDb    = org.Hs.eg.db)
cat("DMR-associated genes mapped to Entrez IDs:", nrow(dmr_gene_entrez), "\n")

# GO enrichment on DEGs
go_deg_bp <- enrichGO(gene          = deg_entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)

go_deg_mf <- enrichGO(gene          = deg_entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)

go_deg_cc <- enrichGO(gene          = deg_entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)

write.csv(as.data.frame(go_deg_bp), "results/enrichment/go_deg_biological_process.csv", row.names = FALSE)
write.csv(as.data.frame(go_deg_mf), "results/enrichment/go_deg_molecular_function.csv", row.names = FALSE)
write.csv(as.data.frame(go_deg_cc), "results/enrichment/go_deg_cellular_component.csv", row.names = FALSE)

dotplot(go_deg_bp, showCategory = 20, title = "GO biological process - DEGs") +
  theme(axis.text.y = element_text(size = 8))
ggsave("results/figures/go_deg_bp_dotplot.png", width = 10, height = 10, dpi = 150)

dotplot(go_deg_mf, showCategory = 20, title = "GO molecular function - DEGs") +
  theme(axis.text.y = element_text(size = 8))
ggsave("results/figures/go_deg_mf_dotplot.png", width = 10, height = 10, dpi = 150)

dotplot(go_deg_cc, showCategory = 20, title = "GO cellular component - DEGs") +
  theme(axis.text.y = element_text(size = 8))
ggsave("results/figures/go_deg_cc_dotplot.png", width = 10, height = 10, dpi = 150)

barplot(go_deg_bp, showCategory = 20, title = "GO biological process - DEGs")
ggsave("results/figures/go_deg_bp_barplot.png", width = 10, height = 10, dpi = 150)

# KEGG pathway enrichment on DEGs
kegg_deg <- enrichKEGG(gene          = deg_entrez$ENTREZID,
                       organism      = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05)

write.csv(as.data.frame(kegg_deg), "results/enrichment/kegg_deg_pathways.csv", row.names = FALSE)

dotplot(kegg_deg, showCategory = 20, title = "KEGG pathways - DEGs") +
  theme(axis.text.y = element_text(size = 8))
ggsave("results/figures/kegg_deg_dotplot.png", width = 10, height = 10, dpi = 150)

# GO and KEGG enrichment on DMR-associated genes
go_dmr_bp <- enrichGO(gene          = dmr_gene_entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)

go_dmr_mf <- enrichGO(gene          = dmr_gene_entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)

write.csv(as.data.frame(go_dmr_bp), "results/enrichment/go_dmr_biological_process.csv", row.names = FALSE)
write.csv(as.data.frame(go_dmr_mf), "results/enrichment/go_dmr_molecular_function.csv", row.names = FALSE)

dotplot(go_dmr_bp, showCategory = 20, title = "GO biological process - DMR genes") +
  theme(axis.text.y = element_text(size = 8))
ggsave("results/figures/go_dmr_bp_dotplot.png", width = 10, height = 10, dpi = 150)

dotplot(go_dmr_mf, showCategory = 20, title = "GO molecular function - DMR genes") +
  theme(axis.text.y = element_text(size = 8))
ggsave("results/figures/go_dmr_mf_dotplot.png", width = 10, height = 10, dpi = 150)

kegg_dmr <- enrichKEGG(gene          = dmr_gene_entrez$ENTREZID,
                       organism      = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05)

write.csv(as.data.frame(kegg_dmr), "results/enrichment/kegg_dmr_pathways.csv", row.names = FALSE)

dotplot(kegg_dmr, showCategory = 20, title = "KEGG pathways - DMR genes") +
  theme(axis.text.y = element_text(size = 8))
ggsave("results/figures/kegg_dmr_dotplot.png", width = 10, height = 10, dpi = 150)

# enrichment map to visualise relationships between GO terms
if (nrow(as.data.frame(go_deg_bp)) > 0) {
  go_deg_bp_sim <- pairwise_termsim(go_deg_bp)
  emapplot(go_deg_bp_sim, showCategory = 30,
           title = "GO BP enrichment map - DEGs")
  ggsave("results/figures/go_deg_bp_emapplot.png", width = 12, height = 10, dpi = 150)
}
