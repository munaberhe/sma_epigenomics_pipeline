#!/usr/bin/env Rscript
# 07i_export_bigwigs.R
# Exports per-condition methylation proportion bigWig files from the pooled cache.
# Muna Berhe -- bt25018 -- QMUL MSc Bioinformatics

.libPaths('/data/home/bt25018/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomeInfoDb)
})
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

METH_CACHE <- 'results/dmr/meth_pooled_cache.rds'
OUT_DIR    <- 'results/bigwigs'
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

CONDITIONS <- c('ASO_VPA', 'ASO_CTRL', 'Scramble_CTRL', 'Scramble_VPA')

# hg38 chromosome lengths (GRCh38.p13)
HG38_SEQLENGTHS <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
  chr21=46709983,  chr22=50818468,  chrX=156040895,  chrY=57227415,
  chrM=16569
)

message('Loading methylation cache...')
meth_pooled <- readRDS(METH_CACHE)
message('Conditions in cache: ', paste(names(meth_pooled), collapse=', '))

for (cond in CONDITIONS) {
  if (!cond %in% names(meth_pooled)) {
    message('SKIP -- not in cache: ', cond)
    next
  }
  message('Exporting: ', cond)
  gr  <- meth_pooled[[cond]]
  cov <- mcols(gr)$readsN
  m   <- mcols(gr)$readsM
  prop <- ifelse(cov >= 4, m / cov, NA_real_)
  keep <- !is.na(prop)
  gr2  <- gr[keep]

  # Keep only standard chromosomes
  std_chrs <- names(HG38_SEQLENGTHS)
  gr2 <- gr2[as.character(seqnames(gr2)) %in% std_chrs]
  seqlevels(gr2) <- std_chrs
  seqlengths(gr2) <- HG38_SEQLENGTHS

  # BigWig needs a score column only
  mcols(gr2) <- NULL
  gr2$score  <- prop[keep][as.character(seqnames(gr[keep])) %in% std_chrs]

  # Sort
  gr2 <- sort(gr2)

  outfile <- file.path(OUT_DIR, paste0(cond, '_methylation.bw'))
  export.bw(gr2, outfile)
  sz <- round(file.size(outfile) / 1e6, 1)
  message('  Saved: ', basename(outfile), ' (', sz, ' MB)')
}
message('Done. BigWigs in: ', OUT_DIR)
