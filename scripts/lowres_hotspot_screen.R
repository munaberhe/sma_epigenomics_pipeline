.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
})
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

CHR_SIZES <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
  chr21=46709983, chr22=50818468, chrX=156040895, chrY=57227415
)

dmr <- readRDS('results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds')
dmr <- dmr[dmr$context == 'CG']
dmr_df <- as.data.frame(dmr)
dmr_df$chr <- as.character(dmr_df$seqnames)
dmr_df <- dmr_df[dmr_df$chr %in% names(CHR_SIZES), ]

BIN <- 50000

stats <- do.call(rbind, lapply(names(CHR_SIZES), function(ch) {
  d <- dmr_df[dmr_df$chr == ch, ]
  chr_len <- CHR_SIZES[ch]
  bins <- seq(0, chr_len, by=BIN)
  bin_counts <- table(cut(d$start, breaks=bins, include.lowest=TRUE))
  data.frame(
    chr            = ch,
    n_dmrs         = nrow(d),
    chr_length_mb  = round(chr_len/1e6, 1),
    dmrs_per_mb    = round(nrow(d) / (chr_len/1e6), 2),
    max_50kb_bin   = max(bin_counts),
    n_hotspot_bins = sum(bin_counts >= 5)
  )
}))

stats <- stats[order(-stats$dmrs_per_mb), ]
message("Top chromosomes by DMR density (ASO_CTRL vs Scramble_CTRL):")
print(stats, row.names=FALSE)
write.csv(stats, 'results/dmr/chr_dmr_density_50kb.csv', row.names=FALSE)
message("Saved: chr_dmr_density_50kb.csv")
