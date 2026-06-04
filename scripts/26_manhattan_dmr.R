#!/usr/bin/env Rscript
.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(ggplot2)
  library(GenomicRanges)
})

setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/manhattan'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

COLS <- c(
  ASO_CTRL      = "#1B4F8A",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#6B7280"
)

CHR_ORDER <- paste0("chr", c(1:22, "X", "Y"))

# Chromosome sizes for cumulative position
chr_sizes <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
  chr21=46709983,  chr22=50818468,  chrX=156040895,  chrY=57227415)

chr_sizes <- chr_sizes[CHR_ORDER]
chr_offsets <- c(0, cumsum(as.numeric(chr_sizes[-length(chr_sizes)])))
names(chr_offsets) <- CHR_ORDER

make_manhattan <- function(dmr_file, contrast_name, cond_colour) {
  dmrs <- readRDS(dmr_file)
  df <- as.data.frame(dmrs)
  df <- df[df$seqnames %in% CHR_ORDER, ]
  df$seqnames <- factor(df$seqnames, levels=CHR_ORDER)
  df$mid <- (df$start + df$end) / 2
  df$cum_pos <- df$mid + chr_offsets[as.character(df$seqnames)]
  df$meth_diff <- mcols(dmrs)$proportion1 - mcols(dmrs)$proportion2
  df$meth_diff <- df$meth_diff[df$seqnames %in% CHR_ORDER]
  df <- df[!is.na(df$cum_pos) & !is.na(df$meth_diff), ]

  # alternate chromosome shading
  df$shade <- as.integer(df$seqnames) %% 2 == 0

  # chromosome midpoints for x-axis labels
  chr_mids <- chr_offsets + chr_sizes/2

  p <- ggplot(df, aes(x=cum_pos, y=meth_diff, colour=shade)) +
    geom_point(size=0.4, alpha=0.6) +
    geom_hline(yintercept=0, colour="grey40", linewidth=0.3) +
    geom_hline(yintercept=c(-0.2, 0.2), linetype="dashed",
               colour="grey60", linewidth=0.3) +
    scale_colour_manual(values=c("FALSE"=cond_colour,
                                 "TRUE"=adjustcolor(cond_colour, alpha.f=0.5)),
                        guide="none") +
    scale_x_continuous(breaks=chr_mids[CHR_ORDER],
                       labels=gsub("chr","",CHR_ORDER),
                       expand=c(0.01, 0)) +
    scale_y_continuous(limits=c(-1, 1),
                       labels=scales::percent_format(1)) +
    labs(title=paste("DMR methylation differences —", contrast_name),
         x="Chromosome", y="Methylation difference (condition 1 - condition 2)") +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold", size=10),
          axis.text.x=element_text(size=7),
          panel.grid.major.y=element_line(colour="grey92"))
  p
}

contrasts <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL", colour=COLS["ASO_CTRL"],
       label="ASO alone vs Scramble CTRL"),
  list(name="Scramble_VPA_vs_Scramble_CTRL", colour=COLS["Scramble_VPA"],
       label="VPA alone vs Scramble CTRL"),
  list(name="ASO_VPA_vs_Scramble_CTRL", colour=COLS["ASO_VPA"],
       label="Combined vs Scramble CTRL"),
  list(name="ASO_VPA_vs_ASO_CTRL", colour=COLS["ASO_VPA"],
       label="VPA effect on ASO background"),
  list(name="ASO_VPA_vs_Scramble_VPA", colour=COLS["ASO_CTRL"],
       label="ASO effect on VPA background")
)

for (ct in contrasts) {
  f <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(f)) { message("missing: ", f); next }
  message("plotting: ", ct$name)
  p <- make_manhattan(f, ct$label, ct$colour)
  out <- file.path(OUT, paste0("manhattan_", ct$name, ".pdf"))
  ggsave(out, p, width=14, height=4, device=cairo_pdf)
  message("  saved: ", basename(out))
}

message("done. outputs in: ", OUT)
