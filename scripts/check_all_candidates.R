.libPaths(c("~/R/library", .libPaths()))
library(GenomicRanges)
library(DMRcaller)
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

# load all four contrasts
aso_alone  <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
vpa_alone  <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
aso_in_vpa <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_VPA.rds")
vpa_in_aso <- readRDS("results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds")

# candidate coords from all_candidate_coords.csv
coords <- read.csv("results/pairwise_context_scan/all_candidate_coords.csv")

cat(sprintf("%-12s %-18s %10s %10s %10s %10s %10s\n",
    "GENE","CATEGORY","ASO_alone","VPA_alone","ASO_in_VPA","VPA_in_ASO","CORRECT"))

for (i in seq_len(nrow(coords))) {
  g   <- coords$name[i]
  cat_label <- coords$category[i]
  gr  <- GRanges(coords$chr[i],
                 IRanges(coords$dmr_start[i], coords$dmr_end[i]))

  n_aso  <- length(subsetByOverlaps(aso_alone,  gr))
  n_vpa  <- length(subsetByOverlaps(vpa_alone,  gr))
  n_aiv  <- length(subsetByOverlaps(aso_in_vpa, gr))
  n_via  <- length(subsetByOverlaps(vpa_in_aso, gr))

  # check category
  correct <- switch(cat_label,
    synergy = n_aso==0 & n_vpa==0 & n_aiv>=1 & n_via>=1,
    ASO_background_specific = n_aso==0 & n_aiv>=1 & n_via==0,
    VPA_background_specific = n_vpa==0 & n_via>=1 & n_aiv==0,
    FALSE
  )

  cat(sprintf("%-12s %-18s %10d %10d %10d %10d %10s\n",
      g, cat_label, n_aso, n_vpa, n_aiv, n_via,
      ifelse(correct, "OK", "*** WRONG ***")))
}
