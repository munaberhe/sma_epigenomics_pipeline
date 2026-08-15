.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

scan_dips <- function(csv_path, contrast_name, min_diff=0.15, min_cpg=4) {
  df <- read.csv(csv_path, stringsAsFactors=FALSE)
  df <- df[df$regionType=="gain",]  # hypo only
  df$meth_diff <- df$proportion2 - df$proportion1
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "" &
           df$meth_diff > min_diff &
           df$cytosinesCount >= min_cpg,]
  df <- df[order(-df$meth_diff),]
  df <- df[!duplicated(df$SYMBOL),]
  df$contrast <- contrast_name
  df
}

aso     <- scan_dips("results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv",
                     "ASO_CTRL")
aso_vpa <- scan_dips("results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv",
                     "ASO_VPA")

cat("\n=== TOP 20 ASO_CTRL dips (ASO alone vs Scramble_CTRL) ===\n")
print(head(aso[,c("seqnames","SYMBOL","GENENAME","meth_diff","annotation")], 20),
      row.names=FALSE)

cat("\n=== TOP 20 ASO_VPA dips (ASO+VPA vs Scramble_CTRL) ===\n")
print(head(aso_vpa[,c("seqnames","SYMBOL","GENENAME","meth_diff","annotation")], 20),
      row.names=FALSE)

# Find genes with big dips in BOTH contrasts - the strongest story
shared <- merge(
  aso[,c("SYMBOL","GENENAME","seqnames","start","end","meth_diff","annotation")],
  aso_vpa[,c("SYMBOL","meth_diff")],
  by="SYMBOL", suffixes=c("_ASO","_ASO_VPA"))
shared <- shared[order(-shared$meth_diff_ASO_VPA),]

cat("\n=== SHARED genes with big dips in BOTH contrasts ===\n")
print(head(shared[,c("seqnames","SYMBOL","GENENAME","meth_diff_ASO","meth_diff_ASO_VPA")], 20),
      row.names=FALSE)

write.csv(shared, "results/dmr_annotation/shared_dip_genes_ASO_and_ASOVPA.csv",
          row.names=FALSE)
message("saved: shared_dip_genes_ASO_and_ASOVPA.csv")
