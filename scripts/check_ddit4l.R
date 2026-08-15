.libPaths(c("~/R/library", .libPaths()))
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

# DDIT4L DMR from annotated file
ann <- read.csv("results/pairwise_context_scan/annotated_ASO_in_VPA.csv")
hits <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL=="DDIT4L",
            c("seqnames","start","end","meth_diff","annotation","pValue")]
cat("DDIT4L DMRs:\n")
print(hits)
