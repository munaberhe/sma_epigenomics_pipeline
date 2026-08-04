#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(UpSetR)
  library(grid)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/figures/upset"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)


hc <- function(gr) gr[mcols(gr)$cytosinesCount >= 6]

aso_alone  <- hc(readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds"))
vpa_alone  <- hc(readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds"))
aso_in_vpa <- hc(readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_VPA.rds"))
vpa_in_aso <- hc(readRDS("results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds"))

sets <- list(
  ASO_alone  = aso_alone,
  VPA_alone  = vpa_alone,
  ASO_in_VPA = aso_in_vpa,
  VPA_in_ASO = vpa_in_aso
)

cat("DMR counts per contrast:\n")
for (n in names(sets)) cat(" ", n, ":", length(sets[[n]]), "\n")


message("Building universe...")
universe <- reduce(do.call(c, unname(sets)))
message("  Universe size: ", length(universe))


message("Computing membership...")
mat <- data.frame(
  ASO_alone  = as.integer(overlapsAny(universe, aso_alone)),
  VPA_alone  = as.integer(overlapsAny(universe, vpa_alone)),
  ASO_in_VPA = as.integer(overlapsAny(universe, aso_in_vpa)),
  VPA_in_ASO = as.integer(overlapsAny(universe, vpa_in_aso))
)

cat("\nMembership column sums (should match DMR counts approximately):\n")
print(colSums(mat))
cat("Universe size:", nrow(mat), "\n")

# Save membership matrix
write.csv(mat, file.path(OUT, "upset_membership_matrix.csv"),
          row.names=FALSE)


message("Plotting UpSet...")
pdf(file.path(OUT, "upset_4contrasts.pdf"),
    width=14, height=7.5, onefile=FALSE)

upset(mat,
  sets           = c("VPA_alone", "VPA_in_ASO", "ASO_in_VPA", "ASO_alone"),
  keep.order     = TRUE,
  order.by       = "freq",
  sets.bar.color = c("#F0A500", "#C0392B", "#B2182B", "#1F3A5F"),
  main.bar.color = "#2C3E50",
  text.scale     = 1.3,
  mb.ratio       = c(0.55, 0.45),
  mainbar.y.label = "DMR regions in this combination",
  sets.x.label    = "DMRs per contrast",
  point.size     = 3,
  line.size      = 1
)

grid.text(
  "DMR overlap across four pairwise contrasts",
  x=0.65, y=0.98,
  gp=gpar(fontsize=12, fontface="bold", col="#1A2A3A")
)
grid.text(
  paste0("Universe = union of all DMRs reduced to non-overlapping intervals",
         " (n=", format(nrow(mat), big.mark=","), ").",
         " cytosinesCount >= 6 filter applied."),
  x=0.65, y=0.945,
  gp=gpar(fontsize=8, col="grey30")
)

dev.off()
message("Saved: upset_4contrasts.pdf")


message("Writing intersection summary...")
get_n <- function(a=0, v=0, aiv=0, via=0) {
  sum(mat$ASO_alone==a & mat$VPA_alone==v &
      mat$ASO_in_VPA==aiv & mat$VPA_in_ASO==via)
}

summary_df <- data.frame(
  intersection = c(
    "ASO alone only",
    "VPA alone only",
    "ASO_in_VPA only",
    "VPA_in_ASO only",
    "ASO alone + VPA alone",
    "ASO alone + ASO_in_VPA",
    "VPA alone + VPA_in_ASO",
    "All four contrasts"
  ),
  n = c(
    get_n(a=1,v=0,aiv=0,via=0),
    get_n(a=0,v=1,aiv=0,via=0),
    get_n(a=0,v=0,aiv=1,via=0),
    get_n(a=0,v=0,aiv=0,via=1),
    get_n(a=1,v=1,aiv=0,via=0),
    get_n(a=1,v=0,aiv=1,via=0),
    get_n(a=0,v=1,aiv=0,via=1),
    get_n(a=1,v=1,aiv=1,via=1)
  )
)

print(summary_df)
write.csv(summary_df, file.path(OUT, "upset_intersection_summary.csv"),
          row.names=FALSE)
message("Done.")
