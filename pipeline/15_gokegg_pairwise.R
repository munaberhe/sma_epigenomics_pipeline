#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/figures/gokegg_pairwise"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

HYPO  <- "#1F3A5F"
HYPER <- "#C0392B"

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       csv="results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv",
       title="ASO alone", short="ASO_alone"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       csv="results/dmr_annotation/Scramble_VPA_vs_Scramble_CTRL_annotated.csv",
       title="VPA alone", short="VPA_alone"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       csv="results/dmr_annotation/ASO_VPA_vs_Scramble_VPA_annotated.csv",
       title="ASO in VPA context", short="ASO_in_VPA"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       csv="results/dmr_annotation/ASO_VPA_vs_ASO_CTRL_annotated.csv",
       title="VPA in ASO context", short="VPA_in_ASO")
)


run_enrichment <- function(ct) {
  message("\n=== ", ct$title, " ===")
  if (!file.exists(ct$csv)) { message("  missing CSV"); return(NULL) }
  df   <- read.csv(ct$csv, stringsAsFactors=FALSE)
  syms <- unique(df$SYMBOL[!is.na(df$SYMBOL) & df$SYMBOL != ""])
  ids  <- bitr(syms, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  if (nrow(ids) < 5) { message("  too few genes"); return(NULL) }

  # GO BP
  go <- tryCatch(
    enrichGO(gene=ids$ENTREZID, OrgDb=org.Hs.eg.db,
             ont="BP", pAdjustMethod="BH",
             pvalueCutoff=0.05, qvalueCutoff=0.2,
             readable=TRUE),
    error=function(e) NULL
  )

  # KEGG
  kk <- tryCatch(
    enrichKEGG(gene=ids$ENTREZID, organism="hsa",
               pvalueCutoff=0.05, pAdjustMethod="BH",
               qvalueCutoff=0.2),
    error=function(e) NULL
  )

  n_go   <- if (!is.null(go)) nrow(as.data.frame(go)) else 0
  n_kegg <- if (!is.null(kk)) nrow(as.data.frame(kk)) else 0
  message("  GO BP terms: ", n_go, " | KEGG terms: ", n_kegg)

  list(go=go, kegg=kk, title=ct$title, short=ct$short, n_genes=nrow(ids))
}

results <- lapply(CONTRASTS, run_enrichment)
names(results) <- sapply(CONTRASTS, function(x) x$short)


make_dotplot <- function(res, type="go", highlight=NULL) {
  enrich <- if (type=="go") res$go else res$kegg
  if (is.null(enrich) || nrow(as.data.frame(enrich)) == 0) {
    p <- ggplot() +
      annotate("text", x=0.5, y=0.5, label="No significant terms",
               size=5, colour="grey50") +
      theme_void() +
      labs(title=paste(res$title, "-", toupper(type)))
    return(p)
  }

  df <- as.data.frame(enrich)

  # highlight pathways that differ between VPA contrasts
  if (!is.null(highlight)) {
    df$highlight <- df$Description %in% highlight
  } else {
    df$highlight <- FALSE
  }

  p <- dotplot(enrich, showCategory=15) +
    ggtitle(paste0(res$title, "\n(", toupper(type), ", n=",
                   res$n_genes, " genes)")) +
    scale_colour_gradient(low=HYPER, high=HYPO,
                          name=expression(-log[10](p.adj))) +
    theme_classic(base_size=10) +
    theme(plot.title     = element_text(face="bold", size=10),
          axis.text.y    = element_text(size=8),
          legend.position = "right")
  p
}

# Save individual PDFs
for (nm in names(results)) {
  r <- results[[nm]]
  if (is.null(r)) next
  # GO
  p_go <- make_dotplot(r, "go")
  ggsave(file.path(OUT, paste0("GO_", nm, ".pdf")),
         p_go, width=9, height=7, device=cairo_pdf)
  # KEGG
  p_kk <- make_dotplot(r, "kegg")
  ggsave(file.path(OUT, paste0("KEGG_", nm, ".pdf")),
         p_kk, width=9, height=7, device=cairo_pdf)
  message("Saved: ", nm)
}


# Find pathways shared between VPA_alone and VPA_in_ASO to highlight
vpa_terms <- if (!is.null(results$VPA_alone$go))
  as.data.frame(results$VPA_alone$go)$Description else character(0)
via_terms  <- if (!is.null(results$VPA_in_ASO$go))
  as.data.frame(results$VPA_in_ASO$go)$Description else character(0)

# Pathways unique to VPA_in_ASO (not in VPA_alone) = ASO modifies VPA effect
vpa_unique <- setdiff(via_terms, vpa_terms)
message("\nPathways unique to VPA_in_ASO (not in VPA_alone): ", length(vpa_unique))
if (length(vpa_unique) > 0) {
  writeLines(vpa_unique, file.path(OUT, "pathways_unique_to_VPA_in_ASO.txt"))
}

go_plots <- lapply(names(results), function(nm) {
  make_dotplot(results[[nm]], "go", highlight=vpa_unique)
})

combined_go <- wrap_plots(go_plots, ncol=2) +
  plot_annotation(
    title   = "GO Biological Process — four pairwise contrasts",
    caption = paste0("Highlighted pathways (if any): unique to VPA-in-ASO context",
                     " (absent from VPA-alone). n=", length(vpa_unique), " such pathways."),
    theme   = theme(
      plot.title   = element_text(face="bold", size=13),
      plot.caption = element_text(colour="grey50", size=9)
    )
  )

ggsave(file.path(OUT, "GO_4contrasts_combined.pdf"),
       combined_go, width=18, height=14, device=cairo_pdf)
ggsave(file.path(OUT, "GO_4contrasts_combined.png"),
       combined_go, width=18, height=14, dpi=150)


kegg_plots <- lapply(names(results), function(nm) {
  make_dotplot(results[[nm]], "kegg")
})

combined_kegg <- wrap_plots(kegg_plots, ncol=2) +
  plot_annotation(
    title = "KEGG pathways — four pairwise contrasts",
    theme = theme(plot.title = element_text(face="bold", size=13))
  )

ggsave(file.path(OUT, "KEGG_4contrasts_combined.pdf"),
       combined_kegg, width=18, height=14, device=cairo_pdf)
ggsave(file.path(OUT, "KEGG_4contrasts_combined.png"),
       combined_kegg, width=18, height=14, dpi=150)

message("\nAll done. Outputs in: ", OUT)
