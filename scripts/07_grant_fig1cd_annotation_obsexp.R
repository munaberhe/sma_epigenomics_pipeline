#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(annotatr)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
})

# stacked-proportion bars + log2(observed/expected) heatmap.
# matches Grant et al 2026 Fig 1C (CpG island context) and Fig 1D (genomic features).
# input : results/dmr/dmr_<contrast>.rds
#         results/dmr/tested_windows_<contrast>.rds  (background = all DMRcaller windows)
# output: results/plots/grant_fig1cd_<contrast>.{pdf,png}

OUT_DIR <- "results/plots"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

CONTRASTS_V <- c(
  "ASO_VPA_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_CTRL_vs_Scramble_CTRL",
  "ASO_VPA_vs_ASO_CTRL",
  "ASO_VPA_vs_Scramble_VPA"
)

# build annotation sets once
cpg_annots <- c(
  "hg38_cpg_islands",
  "hg38_cpg_shores",
  "hg38_cpg_shelves",
  "hg38_cpg_inter"  # open sea
)
feat_annots <- c(
  "hg38_genes_promoters",
  "hg38_genes_5UTRs",
  "hg38_genes_exons",
  "hg38_genes_introns",
  "hg38_genes_3UTRs",
  "hg38_genes_intergenic",
  "hg38_enhancers_fantom"
)

message("building annotations...")
ann_cpg  <- build_annotations(genome="hg38", annotations=cpg_annots)
ann_feat <- build_annotations(genome="hg38", annotations=feat_annots)

# tidy label maps
cpg_label_map <- c(
  "hg38_cpg_islands" = "Island",
  "hg38_cpg_shores"  = "Shore",
  "hg38_cpg_shelves" = "Shelf",
  "hg38_cpg_inter"   = "OpenSea"
)
feat_label_map <- c(
  "hg38_genes_promoters"  = "Promoter",
  "hg38_genes_5UTRs"      = "5'UTR",
  "hg38_genes_exons"      = "Exon",
  "hg38_genes_introns"    = "Intron",
  "hg38_genes_3UTRs"      = "3'UTR",
  "hg38_genes_intergenic" = "Intergenic",
  "hg38_enhancers_fantom" = "Enhancer"
)

# helper: proportion of regions overlapping each annotation type
ann_props <- function(regions, ann, label_map) {
  out <- sapply(names(label_map), function(a) {
    a_gr <- ann[ann$type == a]
    if (length(a_gr) == 0) return(0)
    sum(overlapsAny(regions, a_gr)) / length(regions)
  })
  names(out) <- label_map[names(label_map)]
  out
}

# helper: log2(obs/exp) per annotation
log2_oe <- function(obs_prop, exp_prop) {
  log2((obs_prop + 1e-6) / (exp_prop + 1e-6))
}

build_panel <- function(contrast, ann, label_map, panel_letter, panel_title) {
  rds_path <- file.path("results/dmr", paste0("dmr_", contrast, ".rds"))
  if (!file.exists(rds_path)) return(NULL)
  dmrs <- readRDS(rds_path)
  if (length(dmrs) == 0) return(NULL)

  tested_path <- file.path("results/dmr", paste0("tested_windows_", contrast, ".rds"))
  if (file.exists(tested_path)) {
    bg <- readRDS(tested_path)
  } else {
    # fall back: random sample of equal-size windows across genome
    message("warning: no tested_windows file for ", contrast, ", using DMRs only")
    bg <- dmrs
  }

  obs <- ann_props(dmrs, ann, label_map)
  exp <- ann_props(bg,   ann, label_map)
  oe  <- log2_oe(obs, exp)

  # stacked bar data: DMRs vs background, normalised to 100%
  obs_n <- obs / sum(obs)
  exp_n <- exp / sum(exp)
  bar_df <- bind_rows(
    data.frame(group="Background", category=names(exp_n), prop=exp_n),
    data.frame(group="DMRs",       category=names(obs_n), prop=obs_n)
  )
  bar_df$category <- factor(bar_df$category, levels=label_map)

  # nexus-aligned categorical palette
  pal <- c(
    "Island"="#1B474D","Shore"="#20808D","Shelf"="#BCE2E7","OpenSea"="#28251D",
    "Promoter"="#DA7101","5'UTR"="#FFC553","Exon"="#A84B2F","Intron"="#944454",
    "3'UTR"="#7A39BB","Intergenic"="#848456","Enhancer"="#20808D"
  )

  p_bar <- ggplot(bar_df, aes(x=group, y=prop*100, fill=category)) +
    geom_col(width=0.7, colour="white", linewidth=0.2) +
    scale_fill_manual(values=pal[label_map]) +
    scale_y_continuous(expand=c(0,0), limits=c(0,101)) +
    labs(x=NULL, y="% of regions", fill=NULL,
         title=paste0(panel_letter, "  ", panel_title)) +
    theme_classic(base_size=10) +
    theme(
      legend.position="right",
      legend.key.size=unit(0.35,"cm"),
      legend.text=element_text(size=8),
      plot.title=element_text(size=11, face="bold"),
      axis.text.x=element_text(size=10)
    )

  # log2(obs/exp) heatmap strip
  oe_df <- data.frame(category=names(oe), oe=as.numeric(oe))
  oe_df$category <- factor(oe_df$category, levels=label_map)

  p_oe <- ggplot(oe_df, aes(x=category, y="DMRs", fill=oe)) +
    geom_tile(colour="white") +
    geom_text(aes(label=sprintf("%.2f", oe)), size=3, colour="black") +
    scale_fill_gradient2(low="#006494", mid="#F7F6F2", high="#DA7101",
                         midpoint=0, limits=c(-2,2), oob=scales::squish,
                         name=expression(log[2]~"(obs/exp)")) +
    labs(x=NULL, y=NULL) +
    theme_minimal(base_size=10) +
    theme(
      axis.text.x=element_text(angle=30, hjust=1, size=9),
      panel.grid=element_blank(),
      legend.position="bottom",
      legend.key.height=unit(0.3,"cm"),
      legend.key.width=unit(1,"cm")
    )

  p_bar / p_oe + plot_layout(heights=c(3,1))
}

for (contrast in CONTRASTS_V) {
  message("\nbuilding panels: ", contrast)
  p_c <- build_panel(contrast, ann_cpg,  cpg_label_map,  "C", "CpG island context")
  p_d <- build_panel(contrast, ann_feat, feat_label_map, "D", "Genomic features")
  if (is.null(p_c) || is.null(p_d)) { message("skipping (no data): ", contrast); next }

  combined <- p_c | p_d
  combined <- combined + plot_annotation(
    title = paste0("DMR annotation: ", contrast),
    theme = theme(plot.title = element_text(size=12, face="bold"))
  )

  out_pdf <- file.path(OUT_DIR, paste0("grant_fig1cd_", contrast, ".pdf"))
  out_png <- file.path(OUT_DIR, paste0("grant_fig1cd_", contrast, ".png"))
  ggsave(out_pdf, combined, width=11, height=5)
  ggsave(out_png, combined, width=11, height=5, dpi=300)
  message("wrote: ", out_pdf)
}

message("done.")
