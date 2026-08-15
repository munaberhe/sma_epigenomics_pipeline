#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggplot2)
})

# Annotate DMRs with genomic features and run GO/KEGG enrichment.
# Runs on the genome-wide RDS files produced by 06b_dmrcaller_combine_chr.R.
# Enrichment is split by direction (hypo/hyper) to avoid mixing signals.

OUT_DIR <- "results/dmr_annotation"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

CONTRASTS <- c(
  "ASO_VPA_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_CTRL_vs_Scramble_CTRL",
  "ASO_VPA_vs_ASO_CTRL",
  "ASO_VPA_vs_Scramble_VPA"
)

for (contrast in CONTRASTS) {
  message("\nannotating: ", contrast)
  rds_path <- file.path("results/dmr", paste0("dmr_", contrast, ".rds"))
  if (!file.exists(rds_path)) { message("missing: ", rds_path); next }

  dmrs <- readRDS(rds_path)
  message("DMRs loaded: ", length(dmrs))

  # match chromosome naming to TxDb (UCSC: chr1, chr2...)
  GenomeInfoDb::seqlevelsStyle(dmrs) <- "UCSC"

  # high-confidence filter: at least 6 CpGs per DMR
  dmrs_hc <- dmrs[dmrs$cytosinesCount >= 6]
  message("high-confidence (>=6 CpGs): ", length(dmrs_hc),
          " | hypo: ", sum(dmrs_hc$regionType == "gain"),
          " | hyper: ", sum(dmrs_hc$regionType == "loss"))

  # annotate all DMRs combined
  anno    <- annotatePeak(dmrs_hc, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Hs.eg.db")
  anno_df <- as.data.frame(anno)
  write.csv(anno_df,
            file.path(OUT_DIR, paste0(contrast, "_annotated.csv")),
            row.names=FALSE)

  # top 10 most significant DMRs
  top10 <- anno_df[order(anno_df$pValue), ]
  top10 <- top10[!is.na(top10$SYMBOL), ]
  top10 <- head(top10[, c("seqnames","start","end","cytosinesCount",
                           "regionType","proportion1","proportion2",
                           "pValue","annotation","SYMBOL","GENENAME")], 10)
  top10$direction         <- ifelse(top10$regionType == "gain", "hypo", "hyper")
  top10$methylation_change <- round(top10$proportion2 - top10$proportion1, 3)
  write.csv(top10,
            file.path(OUT_DIR, paste0(contrast, "_top10_genes.csv")),
            row.names=FALSE)
  print(top10[, c("SYMBOL","direction","pValue","methylation_change","annotation")])

  # annotation plots
  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_pie.pdf")), width=8, height=6)
  plotAnnoPie(anno, main=paste0("Genomic features\n", contrast))
  dev.off()

  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_bar.pdf")), width=10, height=6)
  plotAnnoBar(anno, title=paste0("Genomic features — ", contrast))
  dev.off()

  pdf(file.path(OUT_DIR, paste0(contrast, "_TSS_distance.pdf")), width=10, height=6)
  plotDistToTSS(anno, title=paste0("Distance to TSS — ", contrast))
  dev.off()

  # GO and KEGG enrichment — run separately for hypo and hyper
  for (direction in c("gain", "loss")) {
    label <- if (direction == "gain") "hypo" else "hyper"
    sub   <- dmrs_hc[dmrs_hc$regionType == direction]
    message("  ", label, ": ", length(sub), " DMRs")
    if (length(sub) < 10) { message("  too few, skipping"); next }

    anno_sub <- annotatePeak(sub, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
    genes    <- unique(as.data.frame(anno_sub)$geneId)
    genes    <- genes[!is.na(genes)]
    message("  genes: ", length(genes))
    if (length(genes) < 3) { message("  too few genes"); next }

    # GO biological process
    go <- enrichGO(gene=genes, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                   ont="BP", pAdjustMethod="BH",
                   pvalueCutoff=0.05, qvalueCutoff=0.2)
    if (!is.null(go) && nrow(go) > 0) {
      message("  top GO BP: ", go$Description[1])
      write.csv(as.data.frame(go),
                file.path(OUT_DIR, paste0(contrast, "_GO_BP_", label, ".csv")),
                row.names=FALSE)
      pdf(file.path(OUT_DIR, paste0(contrast, "_GO_BP_", label, "_dotplot.pdf")),
          width=10, height=8)
      print(dotplot(go, showCategory=15,
                    title=paste0("GO BP (", label, ") — ", contrast)))
      dev.off()
    } else {
      message("  no significant GO BP terms")
    }

    # KEGG pathways
    kegg <- enrichKEGG(gene=genes, organism="hsa",
                       pvalueCutoff=0.1, qvalueCutoff=0.3, minGSSize=5)
    if (!is.null(kegg) && nrow(kegg) > 0) {
      message("  top KEGG: ", kegg$Description[1])
      write.csv(as.data.frame(kegg),
                file.path(OUT_DIR, paste0(contrast, "_KEGG_", label, ".csv")),
                row.names=FALSE)
      pdf(file.path(OUT_DIR, paste0(contrast, "_KEGG_", label, "_dotplot.pdf")),
          width=10, height=8)
      print(dotplot(kegg, showCategory=15,
                    title=paste0("KEGG (", label, ") — ", contrast)))
      dev.off()
    } else {
      message("  no significant KEGG")
    }
  }
}
message("\ndone. outputs in: ", OUT_DIR)
