#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomeInfoDb)
})

# Extract top 10 hypomethylated DMRs per contrast from annotated CSV files.
# Filters to meth_diff >= 0.20, deduplicates by gene (keeps largest effect),
# then ranks by p-value.

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
OUT  <- "results/dmr_annotation"

CONTRASTS <- c(
  "ASO_CTRL_vs_Scramble_CTRL",
  "ASO_VPA_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_VPA_vs_ASO_CTRL",
  "ASO_VPA_vs_Scramble_VPA"
)

for (contrast in CONTRASTS) {
  message("processing: ", contrast)
  anno_df <- read.csv(file.path(OUT, paste0(contrast, "_annotated.csv")))

  # "gain" = treatment hypomethylated (proportion1 < proportion2)
  hypo           <- anno_df[anno_df$regionType == "gain" & !is.na(anno_df$SYMBOL), ]
  hypo$meth_diff <- hypo$proportion2 - hypo$proportion1
  hypo           <- hypo[hypo$meth_diff >= 0.20, ]

  # deduplicate by gene — keep the DMR with the largest effect per gene
  hypo <- hypo[order(-hypo$meth_diff, hypo$pValue), ]
  hypo <- hypo[!duplicated(hypo$SYMBOL), ]
  hypo <- hypo[order(hypo$pValue), ]

  top10 <- head(hypo[, c("seqnames","start","end","cytosinesCount",
                          "regionType","proportion1","proportion2",
                          "pValue","meth_diff","annotation",
                          "SYMBOL","GENENAME")], 10)
  write.csv(top10,
            file.path(OUT, paste0(contrast, "_top10_hypo.csv")),
            row.names=FALSE)
  print(top10[, c("SYMBOL","meth_diff","pValue","cytosinesCount","annotation")])
}
message("done.")
