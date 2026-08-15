#!/usr/bin/env Rscript
.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/tss_metaplot'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

COLS <- c(ASO_CTRL="#1B4F8A", ASO_VPA="#B2182B",
          Scramble_VPA="#F0A500", Scramble_CTRL="#6B7280")

message("Loading cache...")
cache <- readRDS('results/dmr/meth_pooled_cache.rds')

message("Getting TSS coordinates...")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
genes <- genes[seqnames(genes) %in% paste0("chr", c(1:22,"X"))]
genes <- genes[width(genes) > 1000]

WIN <- 3000
BIN <- 200
n_bins <- WIN*2 / BIN
bin_mids <- seq(-WIN + BIN/2, WIN - BIN/2, by=BIN)

message("Computing TSS profiles (", length(genes), " genes, ", n_bins, " bins)...")

get_tss_profile <- function(meth, genes, win=3000, bin=200) {
  tss_positions <- ifelse(as.character(strand(genes)) == "+",
                          start(genes), end(genes))
  chrs <- as.character(seqnames(genes))
  bins <- seq(-win, win, by=bin)
  
  result <- numeric(length(bins)-1)
  counts <- numeric(length(bins)-1)
  
  meth_df <- data.frame(
    chr = as.character(seqnames(meth)),
    pos = start(meth),
    M   = mcols(meth)$readsM,
    N   = mcols(meth)$readsN
  )
  meth_df <- meth_df[meth_df$N >= 3, ]
  
  for (b in seq_len(length(bins)-1)) {
    bin_start <- bins[b]
    bin_end   <- bins[b+1]
    total_M <- 0
    total_N <- 0
    for (i in seq_along(tss_positions)) {
      tss <- tss_positions[i]
      chr <- chrs[i]
      is_neg <- as.character(strand(genes[i])) == "-"
      sub <- meth_df[meth_df$chr == chr, ]
      if (nrow(sub) == 0) next
      if (is_neg) {
        idx <- sub$pos >= (tss - bin_end) & sub$pos < (tss - bin_start)
      } else {
        idx <- sub$pos >= (tss + bin_start) & sub$pos < (tss + bin_end)
      }
      total_M <- total_M + sum(sub$M[idx])
      total_N <- total_N + sum(sub$N[idx])
    }
    if (total_N > 0) {
      result[b] <- total_M / total_N
      counts[b] <- total_N
    }
  }
  result
}

# Use chr1 only for speed — representative chromosome
message("Subsetting to chr1 for TSS metaplot...")
genes_chr1 <- genes[seqnames(genes) == "chr1"]
message("  genes on chr1: ", length(genes_chr1))

cache_chr1 <- lapply(cache, function(m) m[seqnames(m) == "chr1"])

profiles <- list()
for (cond in names(cache_chr1)) {
  message("  computing: ", cond)
  profiles[[cond]] <- get_tss_profile(cache_chr1[[cond]], genes_chr1,
                                       win=WIN, bin=BIN)
}

prof_df <- do.call(rbind, lapply(names(profiles), function(cond) {
  data.frame(pos=bin_mids, meth=profiles[[cond]], condition=cond)
}))
prof_df$condition <- factor(prof_df$condition,
  levels=c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA"))

p <- ggplot(prof_df, aes(x=pos, y=meth, colour=condition)) +
  geom_line(linewidth=0.9, na.rm=TRUE) +
  geom_vline(xintercept=0, linetype="dashed",
             colour="grey40", linewidth=0.5) +
  annotate("text", x=150, y=max(prof_df$meth, na.rm=TRUE)*0.97,
           label="TSS", size=3.5, colour="grey40", hjust=0) +
  scale_colour_manual(values=COLS) +
  scale_x_continuous(breaks=seq(-3000,3000,1000),
                     labels=function(x) paste0(x/1000,"kb")) +
  scale_y_continuous(labels=scales::percent_format(1)) +
  labs(title="Average CpG methylation around TSS (chr1, protein-coding genes)",
       x="Distance from TSS", y="CpG methylation", colour=NULL) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"),
        legend.position="right")

ggsave(file.path(OUT, "TSS_metaplot.pdf"), p,
       width=9, height=5, device=cairo_pdf)
message("saved: TSS_metaplot.pdf")
message("done.")
