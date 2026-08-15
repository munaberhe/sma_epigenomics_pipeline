.libPaths(c("~/R/library", .libPaths()))
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

mat <- read.csv("results/figures/upset/upset_membership_matrix.csv")

cat("Total union loci:", nrow(mat), "\n\n")

# per-contrast totals in UpSet
cat("Per-contrast totals in UpSet (merged loci):\n")
cat("  ASO alone:", sum(mat$ASO_alone), "\n")
cat("  VPA alone:", sum(mat$VPA_alone), "\n")
cat("  ASO in VPA:", sum(mat$ASO_in_VPA), "\n")
cat("  VPA in ASO:", sum(mat$VPA_in_ASO), "\n\n")

# set-difference sets
aso_context <- mat$ASO_in_VPA == 1 & mat$ASO_alone == 0
vpa_context <- mat$VPA_in_ASO == 1 & mat$VPA_alone == 0

cat("Set-difference derivation:\n")
cat("  ASO_in_VPA loci:", sum(mat$ASO_in_VPA), "\n")
cat("  minus ASO_alone loci also in ASO_in_VPA:",
    sum(mat$ASO_alone == 1 & mat$ASO_in_VPA == 1), "\n")
cat("  = ASO context-dependent:", sum(aso_context), "\n\n")

cat("  VPA_in_ASO loci:", sum(mat$VPA_in_ASO), "\n")
cat("  minus VPA_alone loci also in VPA_in_ASO:",
    sum(mat$VPA_alone == 1 & mat$VPA_in_ASO == 1), "\n")
cat("  = VPA context-dependent:", sum(vpa_context), "\n\n")

cat("Synergy (both context-dependent sets):", sum(aso_context & vpa_context), "\n")
cat("Expected under independence:",
    round(sum(aso_context) * sum(vpa_context) / nrow(mat), 1), "\n")
cat("Enrichment:", round(sum(aso_context & vpa_context) /
    (sum(aso_context) * sum(vpa_context) / nrow(mat)), 2), "x\n")
