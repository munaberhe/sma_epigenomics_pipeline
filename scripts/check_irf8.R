.libPaths(c("~/R/library", .libPaths()))
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
ann <- read.csv("results/pairwise_context_scan/annotated_ASO_in_VPA.csv")
irf8 <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL=="IRF8",
            c("seqnames","start","end","meth_diff","annotation",
              "distanceToTSS","pValue")]
irf8$abs_diff <- abs(irf8$meth_diff)
irf8 <- irf8[order(irf8$abs_diff, decreasing=TRUE),]
cat("IRF8 DMRs (ASO_in_VPA):\n")
print(head(irf8, 5))
cat("\nScoring:\n")
cat("  Promoter annotation: +3\n")
cat("  SMA-relevant list: +3\n")
cat("  |meth_diff| >= 0.40:", irf8$abs_diff[1] >= 0.40, "-> +2\n")
cat("  Total:", 3+3+2, "\n")
