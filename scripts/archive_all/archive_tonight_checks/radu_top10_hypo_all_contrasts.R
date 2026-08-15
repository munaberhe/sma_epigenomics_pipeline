CONTRASTS <- c("ASO_VPA_vs_Scramble_CTRL", "Scramble_VPA_vs_Scramble_CTRL",
               "ASO_CTRL_vs_Scramble_CTRL", "ASO_VPA_vs_ASO_CTRL",
               "ASO_VPA_vs_Scramble_VPA")

for (cn in CONTRASTS) {
  cat("\n==========  ", cn, "  ==========\n")
  csv_path <- file.path("results/dmr_annotation", paste0(cn, "_annotated.csv"))
  if (!file.exists(csv_path)) { cat("  CSV not found\n"); next }
  df <- read.csv(csv_path)
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]

  # hypo only: regionType "gain" means cond1 is lower than cond2
  df <- df[df$regionType == "gain", ]
  df <- df[df$pValue < 0.001, ]
  if (nrow(df) == 0) { cat("  no hits passing filter\n"); next }

  df$methDiff <- df$proportion2 - df$proportion1
  df <- df[order(-df$methDiff), ]
  df_unique <- df[!duplicated(df$SYMBOL), ]

  cat("Top 10 hypomethylated genes, sorted by methylation difference:\n")
  print(head(df_unique[, c("SYMBOL","seqnames","start","methDiff","pValue","annotation")], 10),
        row.names=FALSE)
}
