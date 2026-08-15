.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})
setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/dmr_annotation'

df <- read.csv(paste0(OUT, '/ASO_CTRL_vs_Scramble_CTRL_annotated.csv'))
hypo <- df[df$regionType == 'gain' & !is.na(df$SYMBOL), ]
genes <- unique(hypo$SYMBOL)
eg <- bitr(genes, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)

# KEGG
kegg <- enrichKEGG(gene=eg$ENTREZID, organism='hsa',
                   pvalueCutoff=0.05, qvalueCutoff=0.2)
if (!is.null(kegg) && nrow(kegg) > 0) {
  write.csv(as.data.frame(kegg),
            paste0(OUT, '/ASO_CTRL_vs_Scramble_CTRL_KEGG_hypo.csv'),
            row.names=FALSE)
  p <- dotplot(kegg, showCategory=20) +
    theme_classic(base_size=11) +
    labs(title="KEGG enrichment — ASO hypo DMRs")
  ggsave(paste0(OUT, '/ASO_CTRL_vs_Scramble_CTRL_KEGG_hypo_dotplot.pdf'),
         p, width=10, height=8, device=cairo_pdf)
  message("Saved KEGG hypo dotplot")
} else {
  message("No significant KEGG terms")
}

# Also GO for ASO_VPA vs Scramble_VPA (non-directional, all DMRs)
df2 <- read.csv(paste0(OUT, '/ASO_VPA_vs_Scramble_VPA_annotated.csv'))
all_genes <- unique(df2$SYMBOL[!is.na(df2$SYMBOL)])
eg2 <- bitr(all_genes, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)
go2 <- enrichGO(gene=eg2$ENTREZID, OrgDb=org.Hs.eg.db,
                ont='BP', pvalueCutoff=0.05, qvalueCutoff=0.2,
                readable=TRUE)
if (!is.null(go2) && nrow(go2) > 0) {
  write.csv(as.data.frame(go2),
            paste0(OUT, '/ASO_VPA_vs_Scramble_VPA_GO_BP.csv'),
            row.names=FALSE)
  p2 <- dotplot(go2, showCategory=20) +
    theme_classic(base_size=11) +
    labs(title="GO BP — ASO_VPA vs Scramble_VPA (all DMRs)")
  ggsave(paste0(OUT, '/ASO_VPA_vs_Scramble_VPA_GO_BP_dotplot.pdf'),
         p2, width=10, height=8, device=cairo_pdf)
  message("Saved GO BP dotplot ASO_VPA vs Scramble_VPA")
} else {
  message("No significant GO terms")
}
message("Done")
