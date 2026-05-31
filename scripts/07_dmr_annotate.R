.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggplot2)
})

OUT_DIR <- "results/dmr_annotation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

contrasts <- c("ASO_VPA_vs_Scramble_CTRL",
               "Scramble_VPA_vs_Scramble_CTRL",
               "ASO_CTRL_vs_Scramble_CTRL",
               "ASO_VPA_vs_ASO_CTRL",
               "ASO_VPA_vs_Scramble_VPA")

for (contrast in contrasts) {
  message("\n=== Annotating: ", contrast, " ===")
  rds_path <- file.path("results/dmr", paste0("dmr_", contrast, ".rds"))
  if (!file.exists(rds_path)) { message("Missing: ", rds_path); next }

  dmrs <- readRDS(rds_path)
  message("DMRs loaded: ", length(dmrs))

  # Match chromosome naming to TxDb (UCSC style: chr1, chr2...)
  GenomeInfoDb::seqlevelsStyle(dmrs) <- "UCSC"

  # Filter to high-confidence DMRs only
  dmrs_hc <- dmrs[dmrs$cytosinesCount >= 6]
  message("High-confidence DMRs (>=6 CpGs): ", length(dmrs_hc),
          " | gain (hypo): ", sum(dmrs_hc$regionType == "gain"),
          " | loss (hyper): ", sum(dmrs_hc$regionType == "loss"))

  # --- Combined annotation (all directions) --------------------------------
  anno    <- annotatePeak(dmrs_hc, tssRegion = c(-3000, 3000),
                          TxDb = txdb, annoDb = "org.Hs.eg.db")
  anno_df <- as.data.frame(anno)
  write.csv(anno_df,
            file.path(OUT_DIR, paste0(contrast, "_annotated.csv")),
            row.names = FALSE)

# --- Top 10 most significant DMRs per contrast ---------------------------
  top10 <- anno_df[order(anno_df$pValue), ]
  top10 <- top10[!is.na(top10$SYMBOL), ]
  top10 <- head(top10[, c("seqnames","start","end","cytosinesCount",
                           "regionType","proportion1","proportion2",
                           "pValue","annotation","SYMBOL","GENENAME")], 10)
  top10$direction <- ifelse(top10$regionType == "gain", "hypo", "hyper")
  top10$methylation_change <- round(top10$proportion2 - top10$proportion1, 3)
  write.csv(top10,
            file.path(OUT_DIR, paste0(contrast, "_top10_genes.csv")),
            row.names = FALSE)
  message("Top 10 genes saved for: ", contrast)
  print(top10[, c("SYMBOL","direction","pValue","methylation_change","annotation")])

  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_pie.pdf")), width = 8, height = 6)
  plotAnnoPie(anno, main = paste0("Genomic Feature Distribution\n", contrast))
  dev.off()

  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_bar.pdf")), width = 10, height = 6)
  plotAnnoBar(anno, title = paste0("Genomic Feature Distribution — ", contrast))
  dev.off()

  pdf(file.path(OUT_DIR, paste0(contrast, "_TSS_distance.pdf")), width = 10, height = 6)
  plotDistToTSS(anno, title = paste0("Distance to TSS — ", contrast))
  dev.off()

  message("Saved combined annotation plots")

  # --- Directional GO + KEGG enrichment ------------------------------------
  # Run separately for hypo (gain) and hyper (loss) DMRs.
  # Mixing both directions in one enrichment conflates opposing signals —
  # a gene hypermethylated in one DMR and hypomethylated in another would
  # appear once but represent contradictory biology.
  for (direction in c("gain", "loss")) {
    label <- if (direction == "gain") "hypo" else "hyper"
    sub   <- dmrs_hc[dmrs_hc$regionType == direction]
    message("\n  Direction: ", label, " (", length(sub), " DMRs)")

    if (length(sub) < 10) {
      message("  Too few DMRs — skipping enrichment")
      next
    }

    anno_sub    <- annotatePeak(sub, tssRegion = c(-3000, 3000),
                                TxDb = txdb, annoDb = "org.Hs.eg.db")
    anno_sub_df <- as.data.frame(anno_sub)
    genes       <- unique(anno_sub_df$geneId[!is.na(anno_sub_df$geneId)])
    message("  Genes associated: ", length(genes))

    if (length(genes) < 3) { message("  Too few genes"); next }

    # GO Biological Process
    go <- enrichGO(gene          = genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

    if (!is.null(go) && nrow(go) > 0) {
      message("  Top GO BP: ", go$Description[1], " (", nrow(go), " terms)")
      write.csv(as.data.frame(go),
                file.path(OUT_DIR, paste0(contrast, "_GO_BP_", label, ".csv")),
                row.names = FALSE)
      pdf(file.path(OUT_DIR, paste0(contrast, "_GO_BP_", label, "_dotplot.pdf")),
          width = 10, height = 8)
      print(dotplot(go, showCategory = 15,
                    title = paste0("GO BP (", label, ") — ", contrast)))
      dev.off()
    } else {
      message("  No significant GO BP terms")
    }

    # KEGG pathways
    kegg <- enrichKEGG(gene         = genes,
                       organism     = "hsa",
                       pvalueCutoff = 0.1,
                       qvalueCutoff = 0.3,
                       minGSSize    = 5)

    if (!is.null(kegg) && nrow(kegg) > 0) {
      message("  Top KEGG: ", kegg$Description[1], " (", nrow(kegg), " pathways)")
      write.csv(as.data.frame(kegg),
                file.path(OUT_DIR, paste0(contrast, "_KEGG_", label, ".csv")),
                row.names = FALSE)
      pdf(file.path(OUT_DIR, paste0(contrast, "_KEGG_", label, "_dotplot.pdf")),
          width = 10, height = 8)
      print(dotplot(kegg, showCategory = 15,
                    title = paste0("KEGG (", label, ") — ", contrast)))
      dev.off()
    } else {
      message("  No significant KEGG")
    }
  }
}

message("\nDone. Outputs in: ", OUT_DIR)
