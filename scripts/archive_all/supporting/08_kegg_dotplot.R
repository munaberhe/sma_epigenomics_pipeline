.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(ggplot2); library(patchwork)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/kegg"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

CONTRASTS <- list(
  list(name="ASO_specific",  csv="results/dmr_annotation/ASO_specific_annotated.csv"),
  list(name="VPA",           csv="results/dmr_annotation/Scramble_VPA_vs_Scramble_CTRL_annotated.csv"),
  list(name="ASO_VPA",       csv="results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv")
)

run_kegg <- function(name, csv) {
  message("\n=== ", name, " ===")
  if (!file.exists(csv)) { message("  missing: ", csv); return(NULL) }
  df <- read.csv(csv, stringsAsFactors=FALSE)
  syms <- unique(df$SYMBOL[!is.na(df$SYMBOL) & df$SYMBOL != ""])
  ids <- bitr(syms, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  if (nrow(ids) < 5) { message("  too few mapped"); return(NULL) }
  kk <- enrichKEGG(gene=ids$ENTREZID, organism="hsa",
                   pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
  if (is.null(kk) || nrow(as.data.frame(kk))==0) { message("  no terms"); return(NULL) }
  write.csv(as.data.frame(kk), file.path(OUT_DIR, paste0("kegg_",name,".csv")), row.names=FALSE)
  message("  terms: ", nrow(as.data.frame(kk)))
  list(enrich=kk, contrast=name)
}

render <- function(res) {
  dotplot(res$enrich, showCategory=12) +
    ggtitle(paste0("KEGG pathways: ", res$contrast)) +
    scale_colour_gradient(low="#C0392B", high="#1F3A5F") +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold", size=12),
          axis.text.y=element_text(size=9))
}

results <- list()
for (ct in CONTRASTS) {
  r <- run_kegg(ct$name, ct$csv)
  if (!is.null(r)) {
    results[[ct$name]] <- r
    ggsave(file.path(OUT_DIR, paste0("kegg_dotplot_",ct$name,".pdf")),
           render(r), width=8, height=6, device=cairo_pdf)
    message("saved: kegg_dotplot_", ct$name, ".pdf")
  }
}

if (length(results) >= 2) {
  panels <- lapply(names(results), function(nm) render(results[[nm]]))
  combined <- wrap_plots(panels, ncol=length(panels)) +
    plot_annotation(title="KEGG pathway enrichment at DMR-associated genes",
                    theme=theme(plot.title=element_text(face="bold", size=13)))
  ggsave(file.path(OUT_DIR, "kegg_dotplot_combined.pdf"),
         combined, width=8*length(panels), height=6, device=cairo_pdf)
  message("saved: kegg_dotplot_combined.pdf")
}
message("Done. Outputs in: ", OUT_DIR)
