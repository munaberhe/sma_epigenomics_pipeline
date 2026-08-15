.libPaths(c("~/R/library", .libPaths()))
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
df <- read.csv("results/dmr_benchmark_labelswap_real/labelswap_real_summary.csv")
df <- df[df$method=="bins",]
cat("Label-swap summary (bins method):\n")
print(df[, c("window_size","n_real","n_scr_mean","ratio")])
