SMN2_START <- 70049638
SMN2_END   <- 70078522
FLANK_START <- SMN2_START - 5000
FLANK_END   <- SMN2_END + 25000

for (cn in c("ASO_CTRL_vs_Scramble_CTRL", "Scramble_VPA_vs_Scramble_CTRL",
             "ASO_VPA_vs_Scramble_CTRL")) {
  csv <- file.path("results/dmr_annotation", paste0(cn, "_annotated.csv"))
  if (!file.exists(csv)) {
    cat(sprintf("%-35s: CSV not found\n", cn))
    next
  }
  df <- read.csv(csv)
  df <- df[df$seqnames == "chr5", ]

  gene_body <- df[df$start <= SMN2_END & df$end >= SMN2_START, ]
  flanking  <- df[df$start >= FLANK_START & df$end <= FLANK_END &
                  !(df$start <= SMN2_END & df$end >= SMN2_START), ]

  cat(sprintf("\n=== %s ===\n", cn))
  cat(sprintf("  Gene body DMRs (chr5:%d-%d): %d\n",
              SMN2_START, SMN2_END, nrow(gene_body)))
  if (nrow(gene_body) > 0) print(gene_body[, c("start","end","proportion1","proportion2","pValue")])

  cat(sprintf("  Flanking DMRs (chr5:%d-%d, excluding gene body): %d\n",
              FLANK_START, FLANK_END, nrow(flanking)))
  if (nrow(flanking) > 0) print(flanking[, c("start","end","proportion1","proportion2","pValue")])
}

cat("\n=== Sensitivity scan locus (the +36.3% finding) ===\n")
sens_csv <- "results/smn2_local_dmr/SMN2_sensitive_DMRs_all_contrasts.csv"
if (file.exists(sens_csv)) {
  sens <- read.csv(sens_csv)
  print(sens)
} else {
  cat("File not found at:", sens_csv, "\n")
  cat("Searching for alternative sensitivity scan output...\n")
}
