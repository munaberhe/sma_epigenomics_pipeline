
.libPaths("~/R/library")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

OUT_DIR <- "results/dmr_annotation"

# Reload annotation
rds <- readRDS("results/dmr/ASO_VPA_vs_ASO_CTRL/ASO_VPA_vs_ASO_CTRL_all_chr.rds")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno <- annotatePeak(rds, tssRegion = c(-3000, 3000),
                     TxDb = txdb, annoDb = "org.Hs.eg.db")
anno_df <- as.data.frame(anno)
genes <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
message("Genes: ", length(genes))

# GO with relaxed thresholds
go_bp <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.2,
                  qvalueCutoff  = 0.5)
message("GO BP terms: ", if(is.null(go_bp)) 0 else nrow(go_bp))
if (!is.null(go_bp) && nrow(go_bp) > 0) {
  write.csv(as.data.frame(go_bp),
            file.path(OUT_DIR, "ASO_VPA_vs_ASO_CTRL_GO_BP.csv"), row.names = FALSE)
  pdf(file.path(OUT_DIR, "ASO_VPA_vs_ASO_CTRL_GO_BP_dotplot.pdf"), width = 10, height = 8)
  print(dotplot(go_bp, showCategory = 15,
                title = "GO Biological Process — ASO_VPA vs ASO_CTRL"))
  dev.off()
  message("Saved GO BP dotplot")
}

# KEGG with relaxed thresholds
kegg <- enrichKEGG(gene         = genes,
                   organism     = "hsa",
                   pvalueCutoff = 0.2)
message("KEGG pathways: ", if(is.null(kegg)) 0 else nrow(kegg))
if (!is.null(kegg) && nrow(kegg) > 0) {
  write.csv(as.data.frame(kegg),
            file.path(OUT_DIR, "ASO_VPA_vs_ASO_CTRL_KEGG.csv"), row.names = FALSE)
  pdf(file.path(OUT_DIR, "ASO_VPA_vs_ASO_CTRL_KEGG_dotplot.pdf"), width = 10, height = 8)
  print(dotplot(kegg, showCategory = 15,
                title = "KEGG Pathways — ASO_VPA vs ASO_CTRL"))
  dev.off()
  message("Saved KEGG dotplot")
}
message("Done")
