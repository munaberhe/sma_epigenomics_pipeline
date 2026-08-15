.libPaths(c("~/R/library", .libPaths()))
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
ann <- read.csv("results/pairwise_context_scan/annotated_ASO_in_VPA.csv")
hits <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL=="TMEM179B",
            c("seqnames","start","end","meth_diff","regionType",
              "annotation","distanceToTSS","transcriptId","pValue")]
hits$abs_diff <- abs(hits$meth_diff)
hits <- hits[order(hits$abs_diff, decreasing=TRUE),]
print(hits)
