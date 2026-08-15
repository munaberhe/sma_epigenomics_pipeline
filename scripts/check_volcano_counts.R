.libPaths(c("~/R/library", .libPaths()))
library(GenomicRanges)
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

contrasts <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",     label="ASO alone",    n_table=3423),
  list(name="Scramble_VPA_vs_Scramble_CTRL", label="VPA alone",    n_table=598485),
  list(name="ASO_VPA_vs_Scramble_VPA",       label="ASO in VPA",   n_table=23669),
  list(name="ASO_VPA_vs_ASO_CTRL",           label="VPA in ASO",   n_table=664202)
)

for (ct in contrasts) {
  d <- as.data.frame(readRDS(paste0("results/dmr/dmr_", ct$name, ".rds")))
  d$meth_diff <- d$proportion1 - d$proportion2

  n_total  <- nrow(d)
  n_na     <- sum(is.na(d$meth_diff) | is.na(d$pValue))
  n_exact  <- sum(abs(d$meth_diff) == 0.20, na.rm=TRUE)
  n_hypo   <- sum(d$meth_diff < -0.20 & d$pValue < 0.01, na.rm=TRUE)
  n_hyper  <- sum(d$meth_diff >  0.20 & d$pValue < 0.01, na.rm=TRUE)
  n_sum    <- n_hypo + n_hyper

  cat(sprintf("\n%s:\n", ct$label))
  cat(sprintf("  Table 5.3:    %d\n", ct$n_table))
  cat(sprintf("  RDS total:    %d\n", n_total))
  cat(sprintf("  NAs:          %d\n", n_na))
  cat(sprintf("  Exact 0.20:   %d\n", n_exact))
  cat(sprintf("  Hypo+Hyper:   %d\n", n_sum))
  cat(sprintf("  Deficit:      %d\n", ct$n_table - n_sum))
}
