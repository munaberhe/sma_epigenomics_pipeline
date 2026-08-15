for (cn in c("ASO_VPA_vs_Scramble_CTRL", "ASO_CTRL_vs_Scramble_CTRL",
             "Scramble_VPA_vs_Scramble_CTRL", "ASO_VPA_vs_ASO_CTRL",
             "ASO_VPA_vs_Scramble_VPA")) {
  csv <- file.path("results/dmr_annotation", paste0(cn, "_annotated.csv"))
  if (!file.exists(csv)) {
    cat(sprintf("%-35s: CSV not found\n", cn))
    next
  }
  df <- read.csv(csv)
  hits <- df[!is.na(df$SYMBOL) & df$SYMBOL == "SEMA3C", ]
  if (nrow(hits) > 0) {
    best <- hits[order(hits$pValue), ][1, ]
    cat(sprintf("%-35s: %d DMRs, best p=%.3e, methDiff=%+.3f, regionType=%s\n",
                cn, nrow(hits), best$pValue,
                best$proportion2 - best$proportion1, best$regionType))
  } else {
    cat(sprintf("%-35s: 0 DMRs\n", cn))
  }
}
