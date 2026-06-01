suppressPackageStartupMessages({
  library(UpSetR)
  library(GenomicRanges)
})
.libPaths(c("~/R/library", .libPaths()))
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/dmr_overlap"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

message("Loading DMRs...")
aso     <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
aso_vpa <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
scr_vpa <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")

hc <- function(gr) gr[mcols(gr)$cytosinesCount >= 6]
aso     <- hc(aso)
aso_vpa <- hc(aso_vpa)
scr_vpa <- hc(scr_vpa)

make_membership <- function(sets, names) {
  all_ranges <- reduce(do.call(c, sets))
  mat <- matrix(0L, nrow=length(all_ranges), ncol=length(sets),
                dimnames=list(NULL, names))
  for (i in seq_along(sets)) {
    hits <- findOverlaps(all_ranges, sets[[i]])
    mat[unique(queryHits(hits)), i] <- 1L
  }
  as.data.frame(mat)
}

message("Building membership matrix...")
df <- make_membership(
  list(aso, aso_vpa, scr_vpa),
  c("ASO_CTRL", "ASO_VPA", "Scramble_VPA")
)

# UpSetR needs pdf() opened before calling upset()
# onefile=FALSE prevents the blank first page
pdf(file.path(OUT_DIR, "dmr_upset_plot.pdf"),
    width=12, height=7, onefile=FALSE)
upset(
  df,
  sets           = c("Scramble_VPA", "ASO_VPA", "ASO_CTRL"),
  keep.order     = TRUE,
  order.by       = "freq",
  sets.bar.color = c("#E67E22", "#C0392B", "#2980B9"),
  main.bar.color = "#2C3E50",
  text.scale     = 1.3,
  mb.ratio       = c(0.55, 0.45),
  mainbar.y.label = "DMR intersections (n)",
  sets.x.label    = "DMRs per contrast",
  number.angles  = 0,
  point.size     = 3,
  line.size      = 1
)
grid::grid.text(
  "DMR overlap across three contrasts (cytosinesCount >= 6)",
  x=0.65, y=0.97,
  gp=grid::gpar(fontsize=11, fontface="bold", col="#1A2A3A")
)
dev.off()

message("Saved: dmr_upset_plot.pdf")

overlap_aso_both   <- length(subsetByOverlaps(aso, aso_vpa))
overlap_aso_scrvpa <- length(subsetByOverlaps(aso, scr_vpa))
overlap_both_scr   <- length(subsetByOverlaps(aso_vpa, scr_vpa))
aso_specific       <- length(subsetByOverlaps(
  subsetByOverlaps(aso, aso_vpa), scr_vpa, invert=TRUE))

summary_df <- data.frame(
  comparison = c(
    "ASO_CTRL total (high-conf)",
    "ASO_VPA total (high-conf)",
    "Scramble_VPA total (high-conf)",
    "ASO_CTRL overlapping ASO_VPA",
    "ASO_CTRL overlapping Scramble_VPA",
    "ASO_VPA overlapping Scramble_VPA",
    "ASO-specific (in ASO+ASO_VPA, not Scramble_VPA)"
  ),
  n = c(length(aso), length(aso_vpa), length(scr_vpa),
        overlap_aso_both, overlap_aso_scrvpa, overlap_both_scr,
        aso_specific)
)
print(summary_df, row.names=FALSE)
write.table(summary_df, file.path(OUT_DIR, "dmr_overlap_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
message("Done. Outputs in: ", OUT_DIR)

# 5-contrast upset plot -- supplementary figure showing all contrasts
# The original 3-contrast plot is the main figure; this is for the appendix
message("Building 5-contrast membership matrix...")
aso_ctrl_vpa   <- hc(readRDS("results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds"))
aso_vpa_scrvpa <- hc(readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_VPA.rds"))

df5 <- make_membership(
  list(aso, aso_vpa, scr_vpa, aso_ctrl_vpa, aso_vpa_scrvpa),
  c("ASO_CTRL", "ASO_VPA", "Scramble_VPA", "VPA_on_ASO", "ASO_on_VPA")
)

pdf(file.path(OUT_DIR, "dmr_upset_plot_5contrasts.pdf"),
    width=12, height=8, onefile=FALSE)
upset(
  df5,
  sets       = c("Scramble_VPA","ASO_VPA","ASO_CTRL","VPA_on_ASO","ASO_on_VPA"),
  keep.order = TRUE,
  order.by   = "freq",
  sets.bar.color = c("#E67E22","#C0392B","#2980B9","#8E44AD","#16A085"),
  main.bar.color = "#2C3E50",
  text.scale     = 1.3,
  mb.ratio       = c(0.6, 0.4),
  mainbar.y.label = "DMR intersections (n)",
  sets.x.label    = "DMRs per contrast"
)
grid::grid.text(
  "DMR overlap across all 5 contrasts (high-confidence, cytosinesCount >= 6)",
  x=0.65, y=0.985,
  gp=grid::gpar(fontsize=12, fontface="bold", col="#1A2A3A")
)
grid::grid.text(
  "Supplementary figure -- Muna Berhe, QMUL 2026",
  x=0.65, y=0.965,
  gp=grid::gpar(fontsize=9, col="#6B7C93")
)
dev.off()
message("Saved: dmr_upset_plot_5contrasts.pdf")
