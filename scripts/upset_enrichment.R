.libPaths(c("~/R/library", .libPaths()))
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

mat <- read.csv("results/figures/upset/upset_membership_matrix.csv")

aso_context <- mat$ASO_in_VPA == 1 & mat$ASO_alone == 0
vpa_context <- mat$VPA_in_ASO == 1 & mat$VPA_alone == 0

n_aso  <- as.numeric(sum(aso_context))
n_vpa  <- as.numeric(sum(vpa_context))
n_syn  <- as.numeric(sum(aso_context & vpa_context))
n_union <- as.numeric(nrow(mat))

expected   <- n_aso * n_vpa / n_union
enrichment <- n_syn / expected

cat("ASO context-dependent loci:", n_aso, "\n")
cat("VPA context-dependent loci:", n_vpa, "\n")
cat("Synergy (teal bar):", n_syn, "\n")
cat("Union:", n_union, "\n")
cat("Expected under independence:", round(expected, 1), "\n")
cat("Enrichment:", round(enrichment, 1), "x\n")
