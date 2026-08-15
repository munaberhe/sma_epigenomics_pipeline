genes <- c("CHRNB3", "SLC32A1", "DSC3", "RNF169", "MAP3K15")

for (cn in c("ASO_CTRL_vs_Scramble_CTRL", "Scramble_VPA_vs_Scramble_CTRL",
             "ASO_VPA_vs_Scramble_CTRL")) {
  cat("\n==========  ", cn, "  ==========\n")
  csv_path <- file.path("results/dmr_annotation", paste0(cn, "_annotated.csv"))
  df <- read.csv(csv_path)
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  for (g in genes) {
    hits <- df[df$SYMBOL == g, ]
    if (nrow(hits) == 0) {
      cat(g, ": 0 DMRs\n")
    } else {
      best <- hits[order(hits$pValue), ][1, ]
      cat(sprintf("%-10s: %d DMRs, best p=%.2e, methDiff=%+.3f, regionType=%s\n",
                  g, nrow(hits), best$pValue,
                  best$proportion2 - best$proportion1, best$regionType))
    }
  }
}
