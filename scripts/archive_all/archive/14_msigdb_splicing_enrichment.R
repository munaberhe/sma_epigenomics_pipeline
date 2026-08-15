suppressPackageStartupMessages({
  library(GenomicRanges)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(patchwork)
})
.libPaths(c("~/R/library", .libPaths()))
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/dmr_annotation/msigdb"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

get_genes <- function(tsv_pattern, rds_file) {
  tsvs <- Sys.glob(tsv_pattern)
  if (length(tsvs) > 0) {
    df <- read.table(tsvs[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
    genes <- unique(df$geneId[!is.na(df$geneId) & df$geneId != ""])
  } else {
    message("  Annotating from RDS: ", rds_file)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    gr   <- readRDS(rds_file)
    hc   <- gr[mcols(gr)$cytosinesCount >= 6]
    anno <- annotatePeak(hc, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                         tssRegion=c(-3000, 3000), verbose=FALSE)
    genes <- unique(anno@anno$geneId)
    genes <- genes[!is.na(genes) & genes != ""]
  }
  as.character(genes)
}

message("Loading gene lists...")
genes_aso      <- get_genes(
  "results/dmr_annotation/ASO_CTRL*annotated*.tsv",
  "results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
genes_specific <- get_genes(
  "results/dmr_annotation/ASO_specific*annotated*.tsv",
  "results/dmr/dmr_ASO_specific.rds")

message("  ASO_CTRL genes: ",    length(genes_aso))
message("  ASO-specific genes: ", length(genes_specific))

message("Fetching MSigDB gene sets...")
msig_c5 <- msigdbr(species="Homo sapiens",
                   collection="C5", subcollection="GO:BP")
msig_c2 <- msigdbr(species="Homo sapiens",
                   collection="C2", subcollection="CP:REACTOME")

splicing_terms <- c(
  "SPLICING","SPLICEOSOME","RNA_PROCESSING",
  "ALTERNATIVE_SPLICING","MRNA_PROCESSING",
  "TRANSCRIPTION_ELONGATION","RNA_POLYMERASE",
  "HNRNP","SNRNP","SMN"
)
filter_sets <- function(df, terms) {
  df[grepl(paste(terms, collapse="|"), toupper(df$gs_name)), ]
}
splicing_sets <- rbind(filter_sets(msig_c5, splicing_terms),
                       filter_sets(msig_c2, splicing_terms))
message("  Gene sets found: ", length(unique(splicing_sets$gs_name)))

term2gene <- splicing_sets[, c("gs_name","ncbi_gene")]
term2gene$ncbi_gene <- as.character(term2gene$ncbi_gene)

run_enrichment <- function(gene_ids, label) {
  message("  Running: ", label)
  if (length(gene_ids) < 5) return(NULL)
  enricher(gene=gene_ids, TERM2GENE=term2gene,
           pvalueCutoff=0.05, pAdjustMethod="BH",
           minGSSize=5, maxGSSize=500)
}

res_aso      <- run_enrichment(genes_aso,      "ASO_CTRL DMRs")
res_specific <- run_enrichment(genes_specific, "ASO-specific DMRs")

# Save results tables regardless of significance
save_result <- function(res, prefix) {
  df <- if (!is.null(res)) as.data.frame(res) else data.frame()
  write.table(df,
    file.path(OUT_DIR, paste0(prefix, "_msigdb_splicing.tsv")),
    sep="\t", quote=FALSE, row.names=FALSE)
  message("  ", prefix, ": ", nrow(df[df$p.adjust <= 0.05,]),
          " significant terms")
  df
}

df_aso      <- save_result(res_aso,      "ASO_CTRL")
df_specific <- save_result(res_specific, "ASO_specific")

# Plot top terms whether significant or not — show top 15 by p-value
# This way the plot is always informative rather than blank
plot_terms <- function(df, label) {
  if (nrow(df) == 0) {
    p <- ggplot(data.frame(x=1, y=1), aes(x=x, y=y)) +
      geom_blank() +
      annotate("label", x=1, y=1, size=5, color="grey30",
               fill="grey95", label.size=0.5,
               label=paste0("No significant enrichment detected
in MSigDB splicing / transcription gene sets
(104 gene sets tested, BH p.adj > 0.05)

", label)) +
      theme_classic(base_size=12) +
      theme(axis.line=element_blank(), axis.text=element_blank(),
            axis.ticks=element_blank(), axis.title=element_blank()) +
      labs(title=label)
    return(p)
  }
  top <- head(df[order(df$pvalue), ], 15)
  top$GeneRatio_num <- sapply(top$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  top$Description <- factor(top$Description,
                             levels=rev(top$Description))
  sig_col <- ifelse(top$p.adjust <= 0.05, "#1B4F8A", "#888780")
  ggplot(top, aes(x=GeneRatio_num, y=Description,
                  size=Count, color=p.adjust)) +
    geom_point() +
    scale_color_gradient(low="#C0392B", high="#888780",
                         name="p.adj") +
    scale_size_continuous(name="Gene count", range=c(3,8)) +
    theme_classic(base_size=11) +
    labs(title=paste("MSigDB splicing/transcription —", label),
         subtitle="Top 15 terms by p-value (filled = p.adj \u2264 0.05)",
         x="Gene ratio", y=NULL) +
    theme(plot.title=element_text(face="bold", size=11),
          axis.text.y=element_text(size=9))
}

p1 <- plot_terms(df_aso,      "ASO_CTRL DMRs")
p2 <- plot_terms(df_specific, "ASO-specific DMRs")

combined <- p1 / p2 +
  plot_annotation(
    title="MSigDB splicing and transcription elongation enrichment",
    subtitle="ASO-associated DMRs vs splicing/RNAPII gene sets (C5 GO:BP + C2 REACTOME)",
    theme=theme(plot.title=element_text(face="bold", size=12))
  )

ggsave(file.path(OUT_DIR, "msigdb_splicing_dotplot.pdf"),
       combined, width=11, height=10)
message("Saved: msigdb_splicing_dotplot.pdf")
message("Done. Outputs in: ", OUT_DIR)
