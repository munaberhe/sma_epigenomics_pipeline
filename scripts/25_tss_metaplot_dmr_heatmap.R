#!/usr/bin/env Rscript
.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

COLS <- c(
  ASO_CTRL      = "#1B4F8A",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#6B7280"
)

setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/tss_metaplot'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

message("Loading DMR results and methylation cache...")
cache     <- readRDS('results/dmr/meth_pooled_cache.rds')
dmr_aso   <- readRDS('results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds')

# TSS meta-plot
message("Building TSS meta-plot...")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
genes <- genes[seqnames(genes) %in% paste0("chr", c(1:22, "X"))]

# 4kb window around TSS, 200bp bins
WIN   <- 4000
BIN   <- 200
bins  <- seq(-WIN, WIN, by=BIN)
mids  <- bins[-length(bins)] + BIN/2

get_meta <- function(meth, genes, win=4000, bin=200) {
  bins <- seq(-win, win, by=bin)
  mat  <- matrix(NA, nrow=length(genes), ncol=length(bins)-1)
  for (i in seq_along(genes)) {
    g   <- genes[i]
    tss <- if (as.character(strand(g)) == "+") start(g) else end(g)
    chr <- as.character(seqnames(g))
    sub <- meth[seqnames(meth) == chr &
                start(meth) >= (tss - win) &
                end(meth)   <= (tss + win)]
    if (length(sub) == 0) next
    pos <- start(sub) - tss
    if (as.character(strand(g)) == "-") pos <- -pos
    for (b in seq_len(length(bins)-1)) {
      idx <- pos >= bins[b] & pos < bins[b+1]
      if (sum(idx) > 0) {
        m <- mcols(sub)$sumReadsM[idx]
        n <- mcols(sub)$sumReadsN[idx]
        if (sum(n) > 0) mat[i, b] <- sum(m) / sum(n)
      }
    }
  }
  colMeans(mat, na.rm=TRUE)
}

message("Computing TSS profiles (this takes a few minutes)...")
profiles <- lapply(names(cache), function(cond) {
  message("  ", cond)
  vals <- get_meta(cache[[cond]], genes)
  data.frame(pos=mids, meth=vals, condition=cond)
})
prof_df <- do.call(rbind, profiles)
prof_df$condition <- factor(prof_df$condition,
  levels=c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA"))

p_tss <- ggplot(prof_df, aes(x=pos, y=meth, colour=condition)) +
  geom_line(linewidth=0.9, na.rm=TRUE) +
  geom_vline(xintercept=0, linetype="dashed",
             colour="grey40", linewidth=0.5) +
  annotate("text", x=100, y=max(prof_df$meth, na.rm=TRUE)*0.95,
           label="TSS", size=3.5, colour="grey40", hjust=0) +
  scale_colour_manual(values=COLS) +
  scale_x_continuous(breaks=seq(-4000, 4000, 1000),
                     labels=function(x) paste0(x/1000, "kb")) +
  scale_y_continuous(labels=scales::percent_format(1)) +
  labs(title="Average CpG methylation around TSS",
       x="Distance from TSS", y="CpG methylation",
       colour=NULL) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"),
        legend.position="right")

ggsave(file.path(OUT, "TSS_metaplot.pdf"), p_tss,
       width=9, height=5, device=cairo_pdf)
message("saved: TSS_metaplot.pdf")

# DMR heatmap
message("Building DMR heatmap...")
dmr_top <- dmr_aso[order(mcols(dmr_aso)$pValue)[1:500]]

get_meth_at_dmrs <- function(meth, dmrs) {
  ov <- findOverlaps(dmrs, meth)
  vals <- rep(NA, length(dmrs))
  for (i in seq_along(dmrs)) {
    hits <- subjectHits(ov[queryHits(ov)==i])
    if (length(hits) > 0) {
      m <- sum(mcols(meth)$sumReadsM[hits])
      n <- sum(mcols(meth)$sumReadsN[hits])
      if (n > 0) vals[i] <- m/n
    }
  }
  vals
}

heatmap_mat <- sapply(names(cache), function(cond) {
  message("  extracting methylation for ", cond)
  get_meth_at_dmrs(cache[[cond]], dmr_top)
})
colnames(heatmap_mat) <- names(cache)

# Remove rows with too many NAs
keep <- rowSums(is.na(heatmap_mat)) < 2
heatmap_mat <- heatmap_mat[keep, ]
message("  DMRs with data: ", nrow(heatmap_mat))

col_fun <- colorRamp2(c(0, 0.5, 1),
                      c("#2166AC", "white", "#B2182B"))
col_annot <- HeatmapAnnotation(
  condition = colnames(heatmap_mat),
  col = list(condition = COLS),
  annotation_name_side = "left"
)

pdf(file.path(OUT, "DMR_heatmap_top500.pdf"), width=8, height=12)
Heatmap(heatmap_mat,
        name="CpG meth",
        col=col_fun,
        top_annotation=col_annot,
        show_row_names=FALSE,
        show_column_names=TRUE,
        cluster_rows=TRUE,
        cluster_columns=TRUE,
        column_title="Top 500 ASO DMRs — sample clustering",
        row_title=paste0(nrow(heatmap_mat), " DMRs"),
        use_raster=TRUE)
dev.off()
message("saved: DMR_heatmap_top500.pdf")

message("done. outputs in: ", OUT)
