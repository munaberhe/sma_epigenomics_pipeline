.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")
source("scripts/05_smn2_locus_final.R", local=TRUE)
