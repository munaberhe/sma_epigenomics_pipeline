genes <- c("SMN2", "AGRN", "GRIN2C", "GRIN3A", "NLGN4X", "NLGN2", "SCN7A")

for (cn in c("ASO_VPA_vs_Scramble_CTRL", "ASO_CTRL_vs_Scramble_CTRL",
             "Scramble_VPA_vs_Scramble_CTRL", "ASO_VPA_vs_ASO_CTRL",
             "ASO_VPA_vs_Scramble_VPA")) {
  csv <- file.path("results/dmr_annotation", paste0(cn, "_annotated.csv"))
  if (!file.exists(csv)) next
  df <- read.csv(csv)
  cat("\n==========  ", cn, "  ==========\n")
  for (g in genes) {
    hits <- df[!is.na(df$SYMBOL) & df$SYMBOL == g, ]
    if (nrow(hits) > 0) {
      best <- hits[order(hits$pValue), ][1, ]
      cat(sprintf("%-10s: %d DMRs, p=%.3e, methDiff=%+.3f, %s:%d-%d, %s\n",
                  g, nrow(hits), best$pValue,
                  best$proportion2 - best$proportion1,
                  best$seqnames, best$start, best$end, best$annotation))
    } else {
      cat(sprintf("%-10s: 0 DMRs\n", g))
    }
  }
}
