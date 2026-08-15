.libPaths(c("~/R/library", .libPaths()))
library(GenomicRanges)
library(DMRcaller)
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

aso_alone  <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
vpa_alone  <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
aso_in_vpa <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_VPA.rds")
vpa_in_aso <- readRDS("results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds")

coords <- read.csv("results/pairwise_context_scan/all_candidate_coords.csv")

cat(sprintf("%-12s %10s %10s %10s %10s %20s\n",
    "GENE","ASO_alone","VPA_alone","ASO_in_VPA","VPA_in_ASO","CORRECT_CATEGORY"))

for (i in seq_len(nrow(coords))) {
  g  <- coords$name[i]
  gr <- GRanges(coords$chr[i],
                IRanges(coords$dmr_start[i], coords$dmr_end[i]))

  n_aso <- length(subsetByOverlaps(aso_alone,  gr))
  n_vpa <- length(subsetByOverlaps(vpa_alone,  gr))
  n_aiv <- length(subsetByOverlaps(aso_in_vpa, gr))
  n_via <- length(subsetByOverlaps(vpa_in_aso, gr))

  # derive correct category from data
  correct_cat <- if (n_aso==0 & n_vpa==0 & n_aiv>=1 & n_via>=1) {
    "synergy"
  } else if (n_aso==0 & n_vpa==0 & n_aiv>=1 & n_via==0) {
    "ASO_background_specific"
  } else if (n_vpa==0 & n_aso==0 & n_via>=1 & n_aiv==0) {
    "VPA_background_specific"
  } else {
    paste0("UNCLEAR(a=",n_aso,",v=",n_vpa,",aiv=",n_aiv,",via=",n_via,")")
  }

  cat(sprintf("%-12s %10d %10d %10d %10d %20s\n",
      g, n_aso, n_vpa, n_aiv, n_via, correct_cat))
}
