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

# run GO and KEGG enrichment for a set of gene symbols
enrich_genes <- function(syms, label) {
  syms <- unique(syms[!is.na(syms) & syms != ""])
  if (length(syms) < 5) { message("  too few genes: ", label); return(NULL) }
  ids <- bitr(syms, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  if (nrow(ids) < 5) return(NULL)
  go <- tryCatch(
    enrichGO(gene=ids$ENTREZID, OrgDb=org.Hs.eg.db,
             ont="BP", pAdjustMethod="BH",
             pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE),
    error=function(e) NULL)
  kk <- tryCatch(
    enrichKEGG(gene=ids$ENTREZID, organism="hsa",
               pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2),
    error=function(e) NULL)
  n_go   <- if (!is.null(go)) nrow(as.data.frame(go)) else 0
  n_kegg <- if (!is.null(kk)) nrow(as.data.frame(kk)) else 0
  message("  ", label, " | GO: ", n_go, " | KEGG: ", n_kegg)
  list(go=go, kegg=kk, n_genes=nrow(ids), label=label)
}

run_enrichment <- function(ct) {
  message("\n=== ", ct$title, " ===")
  if (!file.exists(ct$csv)) { message("  missing CSV"); return(NULL) }
  df <- read.csv(ct$csv, stringsAsFactors=FALSE)
  # split by methylation direction
  hypo  <- enrich_genes(df$SYMBOL[df$regionType == "gain"],
                        paste0(ct$title, " (hypomethylated)"))
  hyper <- enrich_genes(df$SYMBOL[df$regionType == "loss"],
                        paste0(ct$title, " (hypermethylated)"))
  list(hypo=hypo, hyper=hyper, title=ct$title, short=ct$short)
}

results <- lapply(CONTRASTS, run_enrichment)
names(results) <- sapply(CONTRASTS, function(x) x$short)

# make a dotplot for one direction result
make_dotplot <- function(dir_res, type="go") {
  empty_plot <- function(label="") {
    ggplot() +
      annotate("text", x=0.5, y=0.5,
               label=ifelse(nchar(label)>0,
                 paste0(label, "\nNo significant terms"),
                 "No significant terms"),
               size=3.5, colour="grey60", hjust=0.5, vjust=0.5) +
      theme_void() +
      theme(plot.margin=margin(2,2,2,2,"mm"))
  }
  if (is.null(dir_res)) return(empty_plot())
  enrich <- if (type=="go") dir_res$go else dir_res$kegg
  if (is.null(enrich) || nrow(as.data.frame(enrich)) == 0) {
    return(NULL)
  }
  # drop degenerate panels with fewer than 5 significant terms
  if (nrow(as.data.frame(enrich)) < 5) {
    message("  Skipping panel with <5 terms: ", dir_res$label)
    return(NULL)
  }
  col <- if (grepl("hypo", dir_res$label)) HYPO else HYPER
  dotplot(enrich, showCategory=12) +
    ggtitle(dir_res$label) +
    scale_colour_gradient(low=adjustcolor(col, 0.5), high=col,
                          name=expression(-log[10](p.adj))) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold", size=11, colour=col),
          axis.text.y=element_text(size=10),
          legend.position="right")
}

# save individual per-contrast PDFs (hypo and hyper side by side)
for (nm in names(results)) {
  r <- results[[nm]]
  if (is.null(r)) next
  p_hypo_go   <- make_dotplot(r$hypo,  "go")
  p_hyper_go  <- make_dotplot(r$hyper, "go")
  p_hypo_kk   <- make_dotplot(r$hypo,  "kegg")
  p_hyper_kk  <- make_dotplot(r$hyper, "kegg")
  go_panels  <- Filter(Negate(is.null), list(p_hypo_go, p_hyper_go))
  kk_panels  <- Filter(Negate(is.null), list(p_hypo_kk, p_hyper_kk))
  if (length(go_panels) == 0) go_panels <- list(ggplot() + theme_void())
  if (length(kk_panels) == 0) kk_panels <- list(ggplot() + theme_void())
  combined_go <- wrap_plots(go_panels, ncol=length(go_panels)) +
    plot_annotation(title=paste0(r$title, " - GO Biological Process"),
                    theme=theme(plot.title=element_text(face="bold", size=11)))
  combined_kk <- wrap_plots(kk_panels, ncol=length(kk_panels)) +
    plot_annotation(title=paste0(r$title, " - KEGG pathways"),
                    theme=theme(plot.title=element_text(face="bold", size=11)))
  ggsave(file.path(OUT, paste0("GO_", nm, ".pdf")),
         combined_go, width=14, height=7, device=cairo_pdf)
  ggsave(file.path(OUT, paste0("KEGG_", nm, ".pdf")),
         combined_kk, width=14, height=7, device=cairo_pdf)
  message("Saved: ", nm)
}

# combined 4-contrast GO panel (8 panels: 4 contrasts x 2 directions)
go_plots <- Filter(Negate(is.null), unlist(lapply(names(results), function(nm) {
  r <- results[[nm]]
  list(make_dotplot(r$hypo, "go"), make_dotplot(r$hyper, "go"))
}), recursive=FALSE))

# determine row heights based on whether panels have content
go_has_content <- sapply(seq(1, length(go_plots), by=2), function(i) {
  left  <- !inherits(go_plots[[i]]$layers[[1]]$geom, "GeomBlank") &
            length(go_plots[[i]]$layers) > 0
  right <- if (i+1 <= length(go_plots))
    length(go_plots[[i+1]]$layers) > 0 else FALSE
  left | right
})
go_heights <- ifelse(go_has_content, 4, 0.6)

go_labels <- LETTERS[seq_along(go_plots)]
go_plots_labelled <- lapply(seq_along(go_plots), function(i) {
  go_plots[[i]] + labs(tag=go_labels[i]) +
    theme(plot.tag=element_text(face="bold", size=14))
})
combined_go <- wrap_plots(go_plots_labelled, ncol=2, heights=rep(go_heights, each=1)) +
  plot_annotation(subtitle="Panels shown where FDR < 0.05 and n >= 5 terms. Direction labelled per panel.",
    theme=theme(plot.title=element_text(face="bold", size=13),
                plot.subtitle=element_text(size=10, colour="grey40"))
  )
ggsave(file.path(OUT, "GO_4contrasts_combined.pdf"),
       combined_go, width=18, height=20, device=cairo_pdf)
ggsave(file.path(OUT, "GO_4contrasts_combined.png"),
       combined_go, width=18, height=20, dpi=150)

# combined 4-contrast KEGG panel
kegg_plots <- Filter(Negate(is.null), unlist(lapply(names(results), function(nm) {
  r <- results[[nm]]
  list(make_dotplot(r$hypo, "kegg"), make_dotplot(r$hyper, "kegg"))
}), recursive=FALSE))

kegg_labels <- LETTERS[seq_along(kegg_plots)]
kegg_plots_labelled <- lapply(seq_along(kegg_plots), function(i) {
  kegg_plots[[i]] + labs(tag=kegg_labels[i]) +
    theme(plot.tag=element_text(face="bold", size=14))
})
combined_kegg <- wrap_plots(kegg_plots_labelled, ncol=2) +
  plot_annotation(subtitle="Panels shown where FDR < 0.05 and n >= 5 terms. Direction labelled per panel.",
    theme=theme(plot.title=element_text(face="bold", size=13),
                plot.subtitle=element_text(size=10, colour="grey40"))
  )
ggsave(file.path(OUT, "KEGG_4contrasts_combined.pdf"),
       combined_kegg, width=18, height=20, device=cairo_pdf)
ggsave(file.path(OUT, "KEGG_4contrasts_combined.png"),
       combined_kegg, width=18, height=20, dpi=150)

# pathways unique to VPA_in_ASO hypomethylated
vpa_terms <- if (!is.null(results$VPA_alone$hypo$go))
  as.data.frame(results$VPA_alone$hypo$go)$Description else character(0)
via_terms  <- if (!is.null(results$VPA_in_ASO$hypo$go))
  as.data.frame(results$VPA_in_ASO$hypo$go)$Description else character(0)
vpa_unique <- setdiff(via_terms, vpa_terms)
message("\nPathways unique to VPA_in_ASO hypo: ", length(vpa_unique))
if (length(vpa_unique) > 0)
  writeLines(vpa_unique, file.path(OUT, "pathways_unique_to_VPA_in_ASO.txt"))

message("\nAll done. Outputs in: ", OUT)
