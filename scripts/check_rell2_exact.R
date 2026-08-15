.libPaths(c("~/R/library", .libPaths()))
library(GenomicRanges)
library(DMRcaller)
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

# exact DMR region
region <- GRanges("chr5", IRanges(141639123, 141639422))

for (ct in c("ASO_CTRL_vs_Scramble_CTRL",
             "Scramble_VPA_vs_Scramble_CTRL",
             "ASO_VPA_vs_Scramble_VPA",
             "ASO_VPA_vs_ASO_CTRL")) {
  dmrs <- readRDS(paste0("results/dmr/dmr_", ct, ".rds"))
  hits <- subsetByOverlaps(dmrs, region)
  cat(ct, ":", length(hits), "DMRs at exact RELL2 bin\n")
}
