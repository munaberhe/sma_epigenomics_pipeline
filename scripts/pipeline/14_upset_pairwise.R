#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(UpSetR)
  library(grid)
  library(ggplot2)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/figures/upset"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# load DMRs — use all DMRs from locked genome-wide parameters
# no additional cytosinesCount filter — already applied during DMR calling
aso_alone  <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
vpa_alone  <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
aso_in_vpa <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_VPA.rds")
vpa_in_aso <- readRDS("results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds")

cat("DMR counts per contrast (should match Table 5.3):\n")
cat("  ASO alone:  ", length(aso_alone),  "\n")
cat("  VPA alone:  ", length(vpa_alone),  "\n")
cat("  ASO in VPA: ", length(aso_in_vpa), "\n")
cat("  VPA in ASO: ", length(vpa_in_aso), "\n")

# Step 1: build DMR-centric universe
# merge all DMRs into non-overlapping loci — one row = one unique genomic locus
message("Building DMR-centric universe (reduce all DMRs to non-overlapping loci)...")
universe <- reduce(c(aso_alone, vpa_alone, aso_in_vpa, vpa_in_aso))
message("  Unique loci: ", length(universe))

# Step 2: membership matrix — one row per locus, one col per contrast
# 1 = this locus is covered by a significant DMR in that contrast
message("Computing locus-level membership...")
mat <- data.frame(
  ASO_alone  = as.integer(overlapsAny(universe, aso_alone)),
  VPA_alone  = as.integer(overlapsAny(universe, vpa_alone)),
  ASO_in_VPA = as.integer(overlapsAny(universe, aso_in_vpa)),
  VPA_in_ASO = as.integer(overlapsAny(universe, vpa_in_aso))
)

# Step 3: sanity check — column sums should approximate DMR counts
# (not exact because reduce() merges overlapping DMRs)
cat("\nColumn sums (loci per contrast — should be close to DMR counts above):\n")
print(colSums(mat))
cat("Total unique loci in union:", nrow(mat), "\n")

write.csv(mat, file.path(OUT, "upset_membership_matrix.csv"), row.names=FALSE)

# Step 4: compute thesis-critical intersections explicitly
get_n <- function(a=0,v=0,aiv=0,via=0) {
  sum(mat$ASO_alone==a & mat$VPA_alone==v &
      mat$ASO_in_VPA==aiv & mat$VPA_in_ASO==via)
}

summary_df <- data.frame(
  intersection = c(
    "ASO alone only",
    "VPA alone only",
    "ASO in VPA only (context-dependent ASO)",
    "VPA in ASO only (context-dependent VPA)",
    "ASO alone + VPA alone (additive)",
    "Synergy: ASO in VPA and VPA in ASO candidates)",
    "VPA alone + VPA in ASO",
    "All four contrasts"
  ),
  n = c(
    get_n(a=1,v=0,aiv=0,via=0),
    get_n(a=0,v=1,aiv=0,via=0),
    get_n(a=0,v=0,aiv=1,via=0),
    get_n(a=0,v=0,aiv=0,via=1),
    get_n(a=1,v=1,aiv=0,via=0),
    get_n(a=0,v=0,aiv=1,via=1),
    get_n(a=0,v=1,aiv=0,via=1),
    get_n(a=1,v=1,aiv=1,via=1)
  )
)
cat("\nThesis-critical intersections:\n")
print(summary_df)
write.csv(summary_df,
          file.path(OUT, "upset_intersection_summary.csv"),
          row.names=FALSE)

# Step 5: plot — highlight thesis-critical intersections
# colour the four thesis-critical bars in red
thesis_intersections <- list(
  c(0,0,1,0),  # ASO in VPA only
  c(0,0,0,1),  # VPA in ASO only
  c(0,0,1,1),  # Synergy: ASO in VPA and VPA in ASO)
  c(1,1,0,0)   # ASO alone + VPA alone (additive)
)

message("Plotting UpSet...")
pdf(file.path(OUT, "upset_4contrasts.pdf"),
    width=14, height=8, onefile=FALSE)
upset(
  mat,
  sets           = c("VPA_in_ASO","VPA_alone","ASO_in_VPA","ASO_alone"),
  keep.order     = TRUE,
  order.by       = "freq",
  nintersects    = 15,
  sets.bar.color = c("#C0392B","#F0A500","#C0392B","#1F3A5F"),
  main.bar.color = "#2C3E50",
  text.scale     = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.1),
  mb.ratio       = c(0.60, 0.40),
  mainbar.y.label = "Unique DMR loci in intersection",
  sets.x.label    = "Total unique DMR loci per contrast",
  point.size     = 3.5,
  line.size      = 1.2,
  show.numbers   = "yes",
  query.legend   = "bottom",
  queries = list(
    list(query=intersects, params=list("ASO_in_VPA"),
         color="#2E86AB", active=TRUE, query.name="ASO in VPA unique"),
    list(query=intersects, params=list("VPA_in_ASO"),
         color="#E84855", active=TRUE, query.name="VPA in ASO unique"),
    list(query=intersects, params=list("ASO_in_VPA","VPA_in_ASO"),
         color="#D19900", active=TRUE, query.name="Pairwise synergy (both combination contrasts)")
  )
)
grid.text(
  "DMR locus overlap across four pairwise contrasts",
  x=0.5, y=0.985,
  gp=gpar(fontsize=13, fontface="bold", col="#1A2A3A")
)
grid.text(
  paste0(
    "n=", format(nrow(mat), big.mark=","),
    " unique loci (GenomicRanges::reduce()). Blue: ASO in VPA unique. ",
    "Red: VPA in ASO unique. Terra: pairwise synergy (both combination contrasts)."
  ),
  x=0.5, y=0.948,
  gp=gpar(fontsize=8, col="grey30")
)
dev.off()
message("Saved: upset_4contrasts.pdf")
message("Done.")
