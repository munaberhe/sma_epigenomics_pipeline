.libPaths(c("~/R/library", .libPaths()))
library(DMRcaller)
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

cache <- readRDS("results/dmr/meth_pooled_cache.rds")
region <- GRanges("chr5", IRanges(70088223, 70088522))

cat("Raw methylation at chr5:70,088,223-70,088,522\n")
cat(sprintf("%-20s %8s %8s %8s %8s\n",
    "Condition","readsM","readsN","prop","nCpGs"))

for (cond in c("Scramble_CTRL","ASO_CTRL","Scramble_VPA","ASO_VPA")) {
  m  <- cache[[cond]]
  m  <- m[m$context == "CG"]
  ov <- subsetByOverlaps(m, region)
  sumM  <- sum(ov$readsM, na.rm=TRUE)
  sumN  <- sum(ov$readsN, na.rm=TRUE)
  prop  <- if (sumN > 0) sumM/sumN else NA
  nCpGs <- length(ov)
  cat(sprintf("%-20s %8d %8d %8.3f %8d\n",
      cond, sumM, sumN, prop, nCpGs))
}
