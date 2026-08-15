genes <- c("SIGMAR1", "AGRN", "GRIN2C", "GRIN3A", "NLGN2", "SCN7A")

# For a real interaction effect, we want to directly compare the raw
# methylation proportions across all four conditions at the same genomic
# window, not just count DMRs. This uses the combo contrast's own
# proportion1/proportion2 columns, which already give us:
#   proportion1 = ASO_VPA (combo)
#   proportion2 = Scramble_CTRL (baseline)
# We need the SAME window's proportions in the OTHER contrasts to see if
# ASO_CTRL and Scramble_VPA sit close to Scramble_CTRL (the additive
# prediction) or whether something more complex is happening.

combo <- read.csv("results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv")
aso   <- read.csv("results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv")
vpa   <- read.csv("results/dmr_annotation/Scramble_VPA_vs_Scramble_CTRL_annotated.csv")

for (g in genes) {
  hit <- combo[!is.na(combo$SYMBOL) & combo$SYMBOL == g, ]
  hit <- hit[order(hit$pValue), ][1, ]
  chr <- hit$seqnames; s <- hit$start; e <- hit$end

  cat(sprintf("\n=== %s (chr%s:%d-%d) ===\n", g, chr, s, e))
  cat(sprintf("  Combo (ASO_VPA):       proportion1=%.3f  (this IS the combo value)\n", hit$proportion1))
  cat(sprintf("  Scramble_CTRL (base):  proportion2=%.3f  (this IS the baseline)\n", hit$proportion2))

  # Find the SAME window in the other two contrasts to see where ASO_CTRL
  # and Scramble_VPA actually sit
  aso_hit <- aso[aso$seqnames == chr & aso$start == s & aso$end == e, ]
  vpa_hit <- vpa[vpa$seqnames == chr & vpa$start == s & vpa$end == e, ]

  if (nrow(aso_hit) > 0) {
    cat(sprintf("  ASO_CTRL at same window:      proportion1=%.3f\n", aso_hit$proportion1[1]))
  } else {
    cat("  ASO_CTRL at same window: NOT in DMR table (no significant call, but data may exist)\n")
  }
  if (nrow(vpa_hit) > 0) {
    cat(sprintf("  Scramble_VPA at same window:  proportion1=%.3f\n", vpa_hit$proportion1[1]))
  } else {
    cat("  Scramble_VPA at same window: NOT in DMR table (no significant call, but data may exist)\n")
  }
}
