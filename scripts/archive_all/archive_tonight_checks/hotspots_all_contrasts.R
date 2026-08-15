CONTRASTS <- c("ASO_VPA_vs_Scramble_CTRL", "Scramble_VPA_vs_Scramble_CTRL",
               "ASO_CTRL_vs_Scramble_CTRL", "ASO_VPA_vs_ASO_CTRL",
               "ASO_VPA_vs_Scramble_VPA")

for (cn in CONTRASTS) {
  cat("\n========== ", cn, " ==========\n")
  csv_path <- file.path("results/dmr_annotation", paste0(cn, "_annotated.csv"))
  if (!file.exists(csv_path)) { cat("  CSV not found\n"); next }
  df <- read.csv(csv_path)
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  if (nrow(df) == 0) { cat("  no annotated DMRs\n"); next }
  df$mb_bin <- paste0(df$seqnames, ":", floor(df$start/1e6))
  tab <- table(df$mb_bin)
  tab <- sort(tab, decreasing=TRUE)
  cat("Top 10 densest 1Mb bins:\n")
  print(head(tab, 10))

  top_bins <- names(head(tab, 5))
  cat("\nGenes in top 5 bins:\n")
  for (b in top_bins) {
    parts <- strsplit(b, ":")[[1]]
    chr <- parts[1]; mb <- as.numeric(parts[2])
    sub <- df[df$seqnames == chr & df$start >= mb*1e6 & df$start < (mb+1)*1e6, ]
    genes <- sort(unique(sub$SYMBOL))
    cat(sprintf("  %s (%d DMRs): %s\n", b, nrow(sub), paste(genes, collapse=", ")))
  }
}
