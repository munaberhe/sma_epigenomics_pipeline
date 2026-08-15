#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggrepel)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/figures/volcano_plots"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

HYPO  <- "#1F3A5F"
HYPER <- "#C0392B"
NS    <- "#7A7974"

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       title="ASO alone\n(ASO_CTRL vs Scramble_CTRL)"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       title="VPA alone\n(Scramble_VPA vs Scramble_CTRL)"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       title="ASO in VPA context\n(ASO_VPA vs Scramble_VPA)"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       title="VPA in ASO context\n(ASO_VPA vs ASO_CTRL)")
)

make_volcano <- function(dmr_df, title, ylim_max) {
  dmr_df <- dmr_df %>%
    mutate(
      meth_diff = proportion1 - proportion2,
      log_p     = -log10(pValue),
      colour    = case_when(
        meth_diff <= -0.20 & pValue < 0.01 ~ "Hypo",
        meth_diff >= 0.20 & pValue < 0.01 ~ "Hyper",
        TRUE ~ "NS"
      )
    )

  n_hypo  <- sum(dmr_df$colour == "Hypo")
  n_hyper <- sum(dmr_df$colour == "Hyper")
  has_hits <- (n_hypo + n_hyper) > 0

  # top 10 by combined rank
  top10 <- dmr_df %>%
    filter(colour != "NS") %>%
    mutate(rank = rank(pValue) * rank(-abs(meth_diff))) %>%
    arrange(rank) %>%
    head(10)

  p <- ggplot(dmr_df, aes(x=meth_diff*100, y=pmin(log_p, ylim_max),
                           colour=colour)) +
    geom_point(size=0.7, alpha=0.6) +
    scale_colour_manual(
      values=c(Hypo=HYPO, Hyper=HYPER, NS=NS),
      labels=c(Hypo=sprintf("Hypo (n=%d)", n_hypo),
               Hyper=sprintf("Hyper (n=%d)", n_hyper),
               NS=sprintf("NS (pre-filtered set; hard limits |diff|=20%%, p=0.01)")),
      name=NULL
    ) +
    geom_vline(xintercept=c(-20, 20), linetype="dashed",
               colour="grey40", linewidth=0.4) +
    geom_hline(yintercept=-log10(0.01), linetype="dashed",
               colour="grey40", linewidth=0.4) +
    coord_cartesian(xlim=c(-100,100), ylim=c(0, ylim_max))

  top3 <- head(top10, 3)
  if (has_hits && nrow(top3) > 0) {
    p <- p + geom_text_repel(
      data=top3,
      aes(x=meth_diff*100, y=pmin(-log10(pValue), ylim_max),
          label=sprintf("%s:%s", seqnames, start)),
      size=3, max.overlaps=10, colour="#1F3B5B", fontface="bold"
    )
  }

  if (!has_hits) {
    p <- p + annotate("text", x=0, y=ylim_max*0.9,
                      label="No hits passing thresholds",
                      size=4, colour="#28251D")
  }

  p <- p +
    labs(title=title,
         x="Methylation difference (%)",
         y=expression(-log[10](p)),
         caption=sprintf("Thresholds: |meth.diff| > 20%%, p < 0.01 | Hypo=%d, Hyper=%d",
                         n_hypo, n_hyper)) +
    theme_classic(base_size=11) +
    theme(
      plot.title      = element_text(face="bold", size=10),
      legend.position = "top",
      panel.grid.major = element_line(colour="grey92", linewidth=0.3),
      plot.caption    = element_text(colour="grey50", size=8)
    )
  p
}

dmr_list <- lapply(CONTRASTS, function(ct) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) return(NULL)
  d <- as.data.frame(readRDS(rds))
  d$seqnames <- as.character(d$seqnames)
  d
})

all_logp <- unlist(lapply(dmr_list, function(d) {
  if (is.null(d)) return(NULL)
  -log10(d$pValue)
}))
ylim_max <- quantile(all_logp, 0.999, na.rm=TRUE)
message("Global y cap: ", round(ylim_max, 2))

plots <- lapply(seq_along(CONTRASTS), function(i) {
  ct <- CONTRASTS[[i]]
  d  <- dmr_list[[i]]
  if (is.null(d)) return(NULL)
  p <- make_volcano(d, ct$title, ylim_max)
  p + labs(tag=c("A","B","C","D")[i]) +
    theme(plot.tag=element_text(face="bold", size=14))
})

combined <- wrap_plots(plots, ncol=2) +
  plot_annotation(caption = "x: methylation difference (treatment - reference). Dashed lines: |meth.diff|=20%, p=0.01. Ceiling at 99.9th percentile of -log10(p); saturation reflects pooled read depth not biological signal.",
    theme   = theme(
      plot.title   = element_text(face="bold", size=13),
      plot.caption = element_text(colour="grey50", size=9)
    )
  )

ggsave(file.path(OUT, "volcano_4contrasts.pdf"),
       combined, width=14, height=12, device=cairo_pdf)
ggsave(file.path(OUT, "volcano_4contrasts.png"),
       combined, width=14, height=12, dpi=150)

for (i in seq_along(CONTRASTS)) {
  ct <- CONTRASTS[[i]]
  if (is.null(plots[[i]])) next
  ggsave(file.path(OUT, paste0("volcano_", ct$name, ".pdf")),
         plots[[i]], width=8, height=6, device=cairo_pdf)
}
message("Done.")
