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

OUT <- "results/figures/significance_plots"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# Palette
HYPO  <- "#1F3A5F"  # terra
HYPER <- "#C0392B"  # gold
NS    <- "#7A7974"  # grey
NAVY  <- "#1F3B5B"

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

CHR_LEVELS <- paste0("chr", c(1:22, "X", "Y"))
CHR_COLS   <- rep(c("#EDECE8", "#D8D7D2"), 12)[1:24]
names(CHR_COLS) <- CHR_LEVELS

# compute cumulative positions for x-axis
chr_sizes <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
  chr21=46709983, chr22=50818468, chrX=156040895, chrY=57227415
)
chr_offsets <- c(0, cumsum(chr_sizes[-length(chr_sizes)]))
names(chr_offsets) <- names(chr_sizes)

make_significance_plot <- function(dmr_df, title, y_cap) {
  dmr_df <- dmr_df %>%
    filter(chr %in% CHR_LEVELS) %>%
    mutate(
      offset   = chr_offsets[chr],
      x_pos    = start + offset,
      log_p    = pmin(-log10(pValue), y_cap),
      colour   = case_when(
        regionType == "loss" ~ HYPO,
        regionType == "gain" ~ HYPER,
        TRUE ~ NS
      ),
      direction_label = case_when(
        regionType == "loss" ~ "Hypo",
        regionType == "gain" ~ "Hyper",
        TRUE ~ "NS"
      )
    )

  # chromosome band background
  chr_bands <- data.frame(
    xmin = chr_offsets[CHR_LEVELS],
    xmax = chr_offsets[CHR_LEVELS] + chr_sizes[CHR_LEVELS],
    fill = CHR_COLS[CHR_LEVELS]
  ) %>% filter(!is.na(xmin))

  # x-axis labels at chromosome midpoints
  chr_mids <- chr_offsets[CHR_LEVELS] + chr_sizes[CHR_LEVELS] / 2
  chr_labels <- gsub("chr", "", CHR_LEVELS)

  # top 5 hits for labelling
  top_hits <- dmr_df %>%
    arrange(pValue) %>%
    head(5)

  p <- ggplot() +
    geom_rect(data=chr_bands,
              aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=fill),
              alpha=0.5, inherit.aes=FALSE) +
    scale_fill_identity() +
    geom_point(data=dmr_df,
               aes(x=x_pos, y=log_p, colour=direction_label),
               size=0.6, alpha=0.7) +
    scale_colour_manual(
      values=c(Hypo=HYPO, Hyper=HYPER, NS=NS),
      labels=c(Hypo="Hypomethylated", Hyper="Hypermethylated", NS="NS"),
      name=NULL
    ) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed",
               colour="grey40", linewidth=0.4) +
    geom_text_repel(data=top_hits,
                    aes(x=x_pos, y=log_p, label=sprintf("%s:%d", chr, start)),
                    size=2.5, max.overlaps=10, colour=NAVY) +
    scale_x_continuous(
      breaks=chr_mids,
      labels=chr_labels,
      expand=expansion(mult=0.01)
    ) +
    coord_cartesian(ylim=c(0, y_cap)) +
    labs(title=title,
         x="Chromosome",
         y=expression(-log[10](p[adj])),
         caption=sprintf("n=%s DMRs | y-axis capped at %.1f (99.9th percentile)",
                         format(nrow(dmr_df), big.mark=","), y_cap)) +
    theme_classic(base_size=11) +
    theme(
      plot.title     = element_text(face="bold", size=10),
      axis.text.x    = element_text(size=7, angle=0),
      legend.position = "top",
      panel.grid.major.y = element_line(colour="grey92", linewidth=0.3),
      plot.caption   = element_text(colour="grey50", size=8)
    )
  p
}

# load all four contrasts
dmr_list <- lapply(CONTRASTS, function(ct) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) return(NULL)
  d <- as.data.frame(readRDS(rds))
  colnames(d)[1] <- "chr"
  d$contrast <- ct$name
  d$title    <- ct$title
  d
})
names(dmr_list) <- sapply(CONTRASTS, function(x) x$name)

# global y cap at 99.9th percentile across all contrasts
all_logp <- unlist(lapply(dmr_list, function(d) {
  if (is.null(d)) return(NULL)
  -log10(d$pValue)
}))
y_cap <- quantile(all_logp, 0.999, na.rm=TRUE)
message("Global y cap: ", round(y_cap, 2))

# make all four panels
plots <- lapply(seq_along(CONTRASTS), function(i) {
  ct <- CONTRASTS[[i]]
  d  <- dmr_list[[ct$name]]
  if (is.null(d)) return(NULL)
  make_significance_plot(d, ct$title, y_cap)
})

# 2x2 combined
combined <- wrap_plots(plots, ncol=2) +
  plot_annotation(
    title   = "Genomic distribution of DMR significance - four pairwise contrasts",
    caption = "Each point = one DMR called by DMRcaller (bins method, strict). Dashed line: p=0.05.",
    theme   = theme(
      plot.title   = element_text(face="bold", size=13),
      plot.caption = element_text(colour="grey50", size=9)
    )
  )

ggsave(file.path(OUT, "dmr_significance_4contrasts.pdf"),
       combined, width=18, height=12, device=cairo_pdf)
ggsave(file.path(OUT, "dmr_significance_4contrasts.png"),
       combined, width=18, height=12, dpi=150)
message("Saved significance plots")

# also save individual panels
for (i in seq_along(CONTRASTS)) {
  ct <- CONTRASTS[[i]]
  if (is.null(plots[[i]])) next
  ggsave(file.path(OUT, paste0("dmr_significance_", ct$name, ".pdf")),
         plots[[i]], width=12, height=5, device=cairo_pdf)
}
message("Done.")
