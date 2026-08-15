#!/usr/bin/env Rscript
# assemble_thesis_figures.R
# Assembles all final thesis figures into numbered PDFs with panel labels.
# Panel label convention: bold uppercase A B C D, top-left, outside plot area.
# Muna Berhe, QMUL 2026

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(pdftools)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(cowplot)
})

SCRATCH <- "/gpfs/scratch/bt25018/sma_epigenomics_pipeline"
OUT     <- file.path(SCRATCH, "results/thesis_figures/final")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# helper: add panel label to a rasterized PDF page
pdf_page_as_grob <- function(pdf_path, page=1) {
  img <- pdftools::pdf_render_page(pdf_path, page=page, dpi=300, numeric=TRUE)
  rasterGrob(img, interpolate=TRUE)
}

label_grob <- function(letter) {
  textGrob(letter, x=0.02, y=0.98, just=c("left","top"),
           gp=gpar(fontface="bold", fontsize=18, fontfamily="sans"))
}

# assemble panel: list of (pdf_path, page, label)
assemble <- function(panels, ncol, output, width=16, height=10) {
  grobs <- lapply(panels, function(p) {
    g <- pdf_page_as_grob(p$path, p$page)
    arrangeGrob(g, label_grob(p$label),
                layout_matrix=matrix(c(2,1), nrow=2),
                heights=unit(c(0.05,0.95), "npc"))
  })
  cairo_pdf(output, width=width, height=height, onefile=TRUE)
  grid.arrange(grobs=grobs, ncol=ncol)
  dev.off()
  message("Saved: ", basename(output))
}

# -------------------------------------------------------------------------
# Fig 5.1 -- QC panel (A=violin B=correlation C=DMR heatmap D=chr1 profile)
assemble(list(
  list(path=file.path(SCRATCH,"results/thesis_figures/Fig5.1b_global_methylation_violin.pdf"),
       page=1, label="A"),
  list(path=file.path(SCRATCH,"results/thesis_figures/Fig5.1c_correlation_heatmap.pdf"),
       page=1, label="B"),
  list(path=file.path(SCRATCH,"results/thesis_figures/Fig_DMR_sample_heatmap.pdf"),
       page=1, label="C"),
  list(path=file.path(SCRATCH,"results/thesis_figures/Fig_lowres_chr1_1Mb.pdf"),
       page=1, label="D")
), ncol=2, output=file.path(OUT,"Fig5.1_QC_panel.pdf"), width=16, height=12)

# Fig 5.2 -- DMR sample heatmap (single panel)
file.copy(file.path(SCRATCH,"results/thesis_figures/Fig_DMR_sample_heatmap.pdf"),
          file.path(OUT,"Fig5.2_DMR_heatmap.pdf"), overwrite=TRUE)
message("Saved: Fig5.2_DMR_heatmap.pdf")

# Fig 5.3 -- diverging bar (single)
file.copy(file.path(SCRATCH,"results/figures/genomic_distribution/chr_dmr_diverging_4contrasts.pdf"),
          file.path(OUT,"Fig5.3_diverging_bar.pdf"), overwrite=TRUE)
message("Saved: Fig5.3_diverging_bar.pdf")

# Fig 5.4 -- volcano (single combined)
file.copy(file.path(SCRATCH,"results/figures/volcano_plots/volcano_4contrasts.pdf"),
          file.path(OUT,"Fig5.4_volcano.pdf"), overwrite=TRUE)
message("Saved: Fig5.4_volcano.pdf")

# Fig 5.5 -- annotation enrichment + metagene (A+B)
assemble(list(
  list(path=file.path(SCRATCH,"results/thesis_figures/Fig_annotation_enrichment_heatmap.pdf"),
       page=1, label="A"),
  list(path=file.path(SCRATCH,"results/thesis_figures/Fig_metagene_profile.pdf"),
       page=1, label="B")
), ncol=2, output=file.path(OUT,"Fig5.5_annotation_metagene.pdf"), width=16, height=7)

# Fig 5.6 -- GO combined (single)
file.copy(file.path(SCRATCH,"results/figures/gokegg_pairwise/GO_4contrasts_combined.pdf"),
          file.path(OUT,"Fig5.6_GO.pdf"), overwrite=TRUE)
message("Saved: Fig5.6_GO.pdf")

# Fig 5.7 -- KEGG combined (single)
file.copy(file.path(SCRATCH,"results/figures/gokegg_pairwise/KEGG_4contrasts_combined.pdf"),
          file.path(OUT,"Fig5.7_KEGG.pdf"), overwrite=TRUE)
message("Saved: Fig5.7_KEGG.pdf")

# Fig 5.8 -- UpSet (single)
file.copy(file.path(SCRATCH,"results/figures/upset/upset_4contrasts.pdf"),
          file.path(OUT,"Fig5.8_UpSet.pdf"), overwrite=TRUE)
message("Saved: Fig5.8_UpSet.pdf")

# Fig 5.9 -- candidates: single relevant contrast per gene (panel 3 = ASO_in_VPA)
# RELL2=synergy(p3), DDIT4L=synergy(p3), TMEM179B=ASO_bg(p3), IRF8=synergy(p3)
assemble(list(
  list(path=file.path(SCRATCH,"results/thesis_figures/locus_candidates/01_synergy/RELL2_4contrasts.pdf"),
       page=3, label="A"),
  list(path=file.path(SCRATCH,"results/thesis_figures/locus_candidates/01_synergy/DDIT4L_4contrasts.pdf"),
       page=3, label="B"),
  list(path=file.path(SCRATCH,"results/thesis_figures/locus_candidates/02_ASO_restricted/TMEM179B_4contrasts.pdf"),
       page=3, label="C"),
  list(path=file.path(SCRATCH,"results/thesis_figures/locus_candidates/01_synergy/IRF8_4contrasts.pdf"),
       page=3, label="D")
), ncol=2, output=file.path(OUT,"Fig5.9_candidates.pdf"), width=16, height=12)

# Fig 5.10 -- SMN2: A=masked four-contrast p4 (VPA_in_ASO) B=extended IGV
assemble(list(
  list(path=file.path(SCRATCH,"results/smn2_locus_final/SMN_locus_masked_ASO_VPA_vs_ASO_CTRL.pdf"),
       page=1, label="A"),
  list(path=file.path(SCRATCH,"results/thesis_figures/smn2_extended_igv/SMN2_extended_IGV_50kb.pdf"),
       page=1, label="B")
), ncol=1, output=file.path(OUT,"Fig5.10_SMN2.pdf"), width=14, height=16)

# -------------------------------------------------------------------------
# Appendix F -- full four-contrast locus plots
for (gene_info in list(
  list(gene="RELL2",    dir="01_synergy",        fignum="F.1"),
  list(gene="DDIT4L",   dir="01_synergy",        fignum="F.2"),
  list(gene="TMEM179B", dir="02_ASO_restricted",  fignum="F.3"),
  list(gene="IRF8",     dir="01_synergy",         fignum="F.4")
)) {
  src <- file.path(SCRATCH, "results/thesis_figures/locus_candidates",
                   gene_info$dir, paste0(gene_info$gene, "_4contrasts.pdf"))
  dst <- file.path(OUT, paste0("Fig", gene_info$fignum, "_",
                                gene_info$gene, "_4contrasts.pdf"))
  file.copy(src, dst, overwrite=TRUE)
  message("Saved: ", basename(dst))
}

message("\nAll thesis figures assembled in: ", OUT)
