#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

OUT_DIR <- "results/qc/for_meeting"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COV_DIR <- "results/alignments/bs/by_chr"

all_groups <- list(
  ASO_VPA       = Sys.glob(file.path(COV_DIR, "ASO_VPA_*_chr13.CpG_report.txt.gz")),
  ASO_CTRL      = Sys.glob(file.path(COV_DIR, "ASO_CTRL_*_chr13.CpG_report.txt.gz")),
  Scramble_VPA  = Sys.glob(file.path(COV_DIR, "Scramble_VPA_*_chr13.CpG_report.txt.gz")),
  Scramble_CTRL = Sys.glob(file.path(COV_DIR, "Scramble_CTRL_*_chr13.CpG_report.txt.gz"))
)

group_cols <- c(
  ASO_VPA       = "#D94F3D",
  ASO_CTRL      = "#2E9B6F",
  Scramble_CTRL = "#1D6FA4",
  Scramble_VPA  = "#F0A500"
)

contrasts <- list(
  list(g1 = "ASO_VPA",      g2 = "Scramble_CTRL", tag = "combined_effect"),
  list(g1 = "Scramble_VPA", g2 = "Scramble_CTRL", tag = "VPA_effect"),
  list(g1 = "ASO_VPA",      g2 = "ASO_CTRL",      tag = "ASO_effect")
)

region      <- GRanges("chr13", IRanges(60e6, 80e6))
window_size <- 50000

message("Reading methylation data for all groups...")
meth_data <- list()
for (nm in names(all_groups)) {
  message("  Pooling: ", nm)
  meth_data[[nm]] <- readBismarkPool(all_groups[[nm]])
}
message("All groups loaded.")

for (ct in contrasts) {
  tag     <- paste0("chr13_hotspot_", ct$tag)
  out_pdf <- file.path(OUT_DIR, paste0("final_lowres_", tag, ".pdf"))
  message("Plotting: ", tag)

  prof1 <- computeMethylationProfile(
    methylationData = meth_data[[ct$g1]],
    region          = region,
    windowSize      = window_size,
    context         = "CG")

  prof2 <- computeMethylationProfile(
    methylationData = meth_data[[ct$g2]],
    region          = region,
    windowSize      = window_size,
    context         = "CG")

  profiles <- GRangesList(
    setNames(list(prof1, prof2), c(ct$g1, ct$g2)))

  pdf(out_pdf, width = 10, height = 4)
  par(mar = c(4, 4, 3, 1) + 0.1)
  plotMethylationProfile(
    methylationProfiles = profiles,
    autoscale  = FALSE,
    title      = sprintf("CpG methylation chr13:60-80Mb (50 kb bins)\n%s vs %s",
                         ct$g1, ct$g2),
    col        = unname(group_cols[c(ct$g1, ct$g2)]),
    pch        = c(15, 16),
    lty        = c(1, 2))
  dev.off()
  message("  Saved: ", out_pdf)
}
message("All done.")
