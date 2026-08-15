.libPaths(c("~/R/library", .libPaths()))
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/dmr_annotation"

CONTRASTS <- c(
  "ASO_CTRL_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_VPA_vs_Scramble_CTRL"
)

for (ct in CONTRASTS) {
  csv <- file.path(OUT_DIR, paste0(ct, "_annotated.csv"))
  if (!file.exists(csv)) { message("missing: ", csv); next }
  
  df <- read.csv(csv, stringsAsFactors=FALSE)
  df$meth_diff <- df$proportion2 - df$proportion1
  
  # hypo only + valid gene symbol
  hypo <- df[df$regionType=="gain" & !is.na(df$SYMBOL) & df$SYMBOL!="",]
  hypo <- hypo[order(-hypo$meth_diff),]
  hypo <- hypo[!duplicated(hypo$SYMBOL),]
  top10 <- head(hypo[, c("SYMBOL","GENENAME","seqnames","start","end",
                          "proportion1","proportion2","meth_diff",
                          "pValue","annotation","cytosinesCount")], 10)
  top10$meth_diff_pct <- round(top10$meth_diff * 100, 1)

  out <- file.path(OUT_DIR, paste0(ct, "_top10_hypo_methdiff.csv"))
  write.csv(top10, out, row.names=FALSE)
  
  cat("\n===", ct, "===\n")
  print(top10[, c("SYMBOL","GENENAME","meth_diff_pct","pValue","annotation")],
        row.names=FALSE)
}
message("\nDone.")

# Bonus: protein-coding named genes only (no LOC, no pseudogenes)
message("\n=== ASO_VPA protein-coding top 10 ===")
df2 <- read.csv("results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv")
df2$meth_diff <- df2$proportion2 - df2$proportion1
hypo2 <- df2[df2$regionType=="gain" &
             !is.na(df2$SYMBOL) & df2$SYMBOL!="" &
             !grepl("^LOC|^LINC|^MIR|antisense|-AS[0-9]", df2$SYMBOL, ignore.case=TRUE) &
             !grepl("pseudogene",
                    df2$SYMBOL, ignore.case=TRUE),]
hypo2 <- hypo2[order(-hypo2$meth_diff),]
hypo2 <- hypo2[!duplicated(hypo2$SYMBOL),]
top10_named <- head(hypo2[,c("SYMBOL","GENENAME","meth_diff","pValue","annotation")], 10)
top10_named$meth_diff_pct <- round(top10_named$meth_diff*100, 1)
print(top10_named[,c("SYMBOL","GENENAME","meth_diff_pct","pValue","annotation")],
      row.names=FALSE)
write.csv(top10_named,
  "results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_top10_hypo_named.csv",
  row.names=FALSE)
