genes <- c("SIGMAR1", "AGRN", "GRIN2C", "GRIN3A", "NLGN4X", "NLGN2", "SCN7A", "SMN2")

cat("=== Annotation quality check: intragenic vs nearest-gene distal calls ===\n\n")
df <- read.csv("results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv")
for (g in genes) {
  hits <- df[!is.na(df$SYMBOL) & df$SYMBOL == g, ]
  if (nrow(hits) == 0) next
  best <- hits[order(hits$pValue), ][1, ]
  is_distal <- grepl("Distal Intergenic", best$annotation)
  cat(sprintf("%-10s: annotation=%-20s distanceToTSS=%-8d %s\n",
              g, best$annotation, best$distanceToTSS,
              if (is_distal) "** NEAREST-GENE CALL, NOT INTRAGENIC **" else "(intragenic/promoter)"))
}

cat("\n=== Sanity check: is the filter itself biased? ===\n")
cat("Testing how many genes pass the filter PURELY BY CHANCE\n")
cat("by checking the SAME logic against a clearly-uninteresting random gene set\n\n")

# Pull total gene counts from each list to sanity check the filter denominator
combo <- read.csv("results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv")
aso   <- read.csv("results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv")
vpa   <- read.csv("results/dmr_annotation/Scramble_VPA_vs_Scramble_CTRL_annotated.csv")

combo_genes <- unique(combo$SYMBOL[!is.na(combo$SYMBOL) & combo$SYMBOL != ""])
aso_genes   <- unique(aso$SYMBOL[!is.na(aso$SYMBOL) & aso$SYMBOL != ""])
vpa_genes   <- unique(vpa$SYMBOL[!is.na(vpa$SYMBOL) & vpa$SYMBOL != ""])

cat("Total unique genes in combination contrast:", length(combo_genes), "\n")
cat("Total unique genes in ASO-alone contrast:  ", length(aso_genes), "\n")
cat("Total unique genes in VPA-alone contrast:  ", length(vpa_genes), "\n")
cat("Genes in combo but NOT in ASO-alone:        ", sum(!combo_genes %in% aso_genes), "\n")
cat("Genes in combo but NOT in VPA-alone:        ", sum(!combo_genes %in% vpa_genes), "\n")
cat("Genes in combo but NOT in EITHER:           ",
    sum(!combo_genes %in% aso_genes & !combo_genes %in% vpa_genes), "\n")
