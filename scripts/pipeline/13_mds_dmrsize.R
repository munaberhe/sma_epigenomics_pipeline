#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/figures/qc_plots"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

COND_COLS <- c(
  ASO_CTRL      = "#1B4F8A",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500"
)


message("Building MDS plot...")
cor_csv <- "results/qc/methylation/methylation_correlation.csv"
if (file.exists(cor_csv)) {
  cor_mat <- as.matrix(read.csv(cor_csv, row.names=1))
  # convert correlation to distance
  dist_mat <- as.dist(1 - cor_mat)
  mds      <- cmdscale(dist_mat, k=2)
  mds_df   <- data.frame(
    sample    = rownames(mds),
    MDS1      = mds[,1],
    MDS2      = mds[,2]
  )
  # extract condition from sample name
  mds_df$condition <- gsub("_[123]$", "", mds_df$sample)
  mds_df$replicate <- gsub(".*_([123])$", "\\1", mds_df$sample)
  mds_df$condition <- factor(mds_df$condition,
                             levels=names(COND_COLS))

  p_mds <- ggplot(mds_df, aes(x=MDS1, y=MDS2,
                               colour=condition, label=sample)) +
    geom_point(size=4, alpha=0.9) +
    geom_text(nudge_y=0.005, size=2.8, show.legend=FALSE) +
    scale_colour_manual(values=COND_COLS, name="Condition") +
    stat_ellipse(aes(group=condition), level=0.8,
                 linetype="dashed", linewidth=0.5) +
    labs(title="MDS plot - CpG methylation profiles",
         subtitle="Distance = 1 - Pearson correlation of genome-wide methylation",
         x="MDS dimension 1", y="MDS dimension 2") +
    theme_classic(base_size=12) +
    theme(plot.title    = element_text(face="bold"),
          legend.position = "right",
          panel.grid.major = element_line(colour="grey92", linewidth=0.3))

  ggsave(file.path(OUT, "mds_methylation.pdf"),
         p_mds, width=8, height=6, device=cairo_pdf)
  ggsave(file.path(OUT, "mds_methylation.png"),
         p_mds, width=8, height=6, dpi=150)
  message("Saved MDS plot")
} else {
  message("WARNING: correlation CSV not found at ", cor_csv)
}


message("Building DMR size distribution...")

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",    title="ASO alone"),
  list(name="Scramble_VPA_vs_Scramble_CTRL", title="VPA alone"),
  list(name="ASO_VPA_vs_Scramble_VPA",       title="ASO in VPA context"),
  list(name="ASO_VPA_vs_ASO_CTRL",           title="VPA in ASO context")
)

HYPO  <- "#1F3A5F"
HYPER <- "#C0392B"

size_plots <- lapply(CONTRASTS, function(ct) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) return(NULL)
  dmrs <- as.data.frame(readRDS(rds))
  dmrs$direction <- ifelse(dmrs$regionType == "loss", "Hypo", "Hyper")
  dmrs$width_kb  <- dmrs$width / 1000

  # cap at 10kb for readability
  dmrs_plot <- dmrs[dmrs$width_kb <= 10, ]
  pct_capped <- round(100 * (1 - nrow(dmrs_plot)/nrow(dmrs)), 1)

  ggplot(dmrs_plot, aes(x=width_kb, fill=direction)) +
    geom_histogram(bins=50, colour="white", linewidth=0.1, alpha=0.85) +
    scale_fill_manual(values=c(Hypo=HYPO, Hyper=HYPER)) +
    scale_x_continuous(breaks=seq(0,10,2)) +
    labs(title=ct$title,
         x="DMR width (kb)", y="Count", fill=NULL,
         caption=sprintf("n=%s DMRs | %.1f%% >10kb not shown",
                         format(nrow(dmrs), big.mark=","), pct_capped)) +
    theme_classic(base_size=12) +
    theme(plot.title    = element_text(face="bold", size=10),
          legend.position = "top",
          plot.caption  = element_text(colour="grey50", size=8),
          panel.grid.major.y = element_line(colour="grey92", linewidth=0.3))
})

combined_size <- wrap_plots(Filter(Negate(is.null), size_plots), ncol=2) +
  plot_annotation(
    title   = "DMR size distribution - four pairwise contrasts",
    caption = "Each DMR = one 300bp window passing significance threshold.",
    theme   = theme(plot.title = element_text(face="bold", size=13))
  )

ggsave(file.path(OUT, "dmr_size_distribution.pdf"),
       combined_size, width=12, height=8, device=cairo_pdf)
ggsave(file.path(OUT, "dmr_size_distribution.png"),
       combined_size, width=12, height=8, dpi=150)
message("Saved DMR size distribution")
message("All done.")
