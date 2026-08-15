
.libPaths("~/R/library")
# dmr_annotate.R
# DMR Annotation, Pathway Enrichment and Visualisation
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DMRcaller)
library(ggplot2)

OUT_DIR <- "results/dmr_annotation_0p3"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

contrasts <- c("ASO_VPA_vs_ASO_CTRL", "ASO_VPA_vs_Scramble_CTRL", "VPA_vs_Scramble_CTRL")

for (contrast in contrasts) {
  message("\nAnnotating: ", contrast)
  rds_path <- file.path("results/dmr_0p3", contrast, paste0(contrast, "_all_chr.rds"))
  if (!file.exists(rds_path)) { message("Missing: ", rds_path); next }

  dmrs <- readRDS(rds_path)
  message("  DMRs loaded: ", length(dmrs))

  anno <- annotatePeak(dmrs, tssRegion = c(-3000, 3000),
                       TxDb = txdb, annoDb = "org.Hs.eg.db")

  anno_df <- as.data.frame(anno)
  write.csv(anno_df, file.path(OUT_DIR, paste0(contrast, "_annotated.csv")), row.names = FALSE)
  message("  Saved annotation CSV")

  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_pie.pdf")), width = 8, height = 6)
  plotAnnoPie(anno, main = paste0("Genomic Feature Distribution\n", contrast))
  dev.off()
  message("  Saved pie chart")

  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_bar.pdf")), width = 10, height = 6)
  plotAnnoBar(anno, title = paste0("Genomic Feature Distribution — ", contrast))
  dev.off()
  message("  Saved bar chart")

  pdf(file.path(OUT_DIR, paste0(contrast, "_TSS_distance.pdf")), width = 10, height = 6)
  plotDistToTSS(anno, title = paste0("Distance to TSS — ", contrast))
  dev.off()
  message("  Saved TSS distance plot")

  genes <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
  message("  Genes associated with DMRs: ", length(genes))

  if (length(genes) >= 3) {
    go_bp <- enrichGO(gene          = genes,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)
    if (!is.null(go_bp) && nrow(go_bp) > 0) {
      write.csv(as.data.frame(go_bp),
                file.path(OUT_DIR, paste0(contrast, "_GO_BP.csv")), row.names = FALSE)
      pdf(file.path(OUT_DIR, paste0(contrast, "_GO_BP_dotplot.pdf")), width = 10, height = 8)
      print(dotplot(go_bp, showCategory = 15,
                    title = paste0("GO Biological Process — ", contrast)))
      dev.off()
      message("  GO BP: ", nrow(go_bp), " terms found")
    } else {
      message("  No significant GO BP terms")
    }

    kegg <- enrichKEGG(gene         = genes,
                       organism     = "hsa",
                       pvalueCutoff = 0.05)
    if (!is.null(kegg) && nrow(kegg) > 0) {
      write.csv(as.data.frame(kegg),
                file.path(OUT_DIR, paste0(contrast, "_KEGG.csv")), row.names = FALSE)
      pdf(file.path(OUT_DIR, paste0(contrast, "_KEGG_dotplot.pdf")), width = 10, height = 8)
      print(dotplot(kegg, showCategory = 15,
                    title = paste0("KEGG Pathways — ", contrast)))
      dev.off()
      message("  KEGG: ", nrow(kegg), " pathways found")
    } else {
      message("  No significant KEGG pathways")
    }
  } else {
    message("  Too few genes for enrichment (n=", length(genes), ")")
  }
}

message("\nDone. Outputs in: ", OUT_DIR)
