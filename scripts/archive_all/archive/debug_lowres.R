.libPaths('~/R/library')
library(DMRcaller)
library(GenomicRanges)
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

cache <- readRDS('results/dmr/meth_pooled_cache.rds')
region_chr5 <- GRanges("chr5", IRanges(69000000, 71500000))
region_chr1 <- GRanges("chr1", IRanges(1, 248956422))

prof_chr5 <- computeMethylationProfile(cache[["ASO_CTRL"]], region_chr5,
                                        windowSize=10000, context="CG")
prof_chr1 <- computeMethylationProfile(cache[["ASO_CTRL"]], region_chr1,
                                        windowSize=1000000, context="CG")

cat("chr5 region bins:", length(prof_chr5), "\n")
cat("chr5 positions range:", min(start(prof_chr5)), "-", max(end(prof_chr5)), "\n")
cat("chr1 bins:", length(prof_chr1), "\n")
cat("chr1 positions range:", min(start(prof_chr1)), "-", max(end(prof_chr1)), "\n")

cat("chr5 mean methylation:", round(mean(prof_chr5$sumReadsM/prof_chr5$sumReadsN, na.rm=TRUE), 3), "\n")
cat("chr1 mean methylation:", round(mean(prof_chr1$sumReadsM/prof_chr1$sumReadsN, na.rm=TRUE), 3), "\n")
