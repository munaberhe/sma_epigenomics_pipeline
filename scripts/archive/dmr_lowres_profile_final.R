#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

OUT_DIR <- "results/qc/for_meeting"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BS_DIR <- "results/alignments/bs"

# Group file lists
group_files <- list(
  ASO_VPA       = Sys.glob(file.path(BS_DIR, "ASO_VPA_*_CX_report.txt.CpG_report.txt.gz")),
  ASO_CTRL      = Sys.glob(file.path(BS_DIR, "ASO_CTRL_*_CX_report.txt.CpG_report.txt.gz")),
  Scramble_VPA  = Sys.glob(file.path(BS_DIR, "Scramble_VPA_*_CX_report.txt.CpG_report.txt.gz")),
  Scramble_CTRL = Sys.glob(file.path(BS_DIR, "Scramble_CTRL_*_CX_report.txt.CpG_report.txt.gz"))
)

# Two regions: full chr1 (1Mb bins) and DMR-dense zoom (50kb bins)
regions_list <- list(
  list(region = GRanges("chr1", IRanges(1, 2.5e8)),
       window = 1000000, label = "chr1_1Mb"),
  list(region = GRanges("chr1", IRanges(210e6, 230e6)),
       window = 50000,   label = "chr1_210_230Mb_50kb")
)

# Contrasts to plot
contrasts <- list(
  list(g1 = "ASO_VPA",      g2 = "Scramble_CTRL", tag = "combined_effect"),
  list(g1 = "Scramble_VPA", g2 = "Scramble_CTRL", tag = "VPA_effect")
)

# Group colours — consistent palette
# One colour per group — consistent across all plots in the deck
group_cols <- c(
  ASO_VPA       = "#D94F3D",   # brick red
  ASO_CTRL      = "#2E9B6F",   # teal green
  Scramble_CTRL = "#1D6FA4",   # steel blue
  Scramble_VPA  = "#F0A500"    # amber
)

# Read and pool each group (do this once)
message("Reading methylation data for all groups...")
meth_data <- list()
for (nm in names(group_files)) {
  message("  Pooling: ", nm)
  meth_data[[nm]] <- readBismarkPool(group_files[[nm]])
}
message("All groups loaded.")

# Generate plots
for (ct in contrasts) {
  for (reg in regions_list) {
    tag     <- paste0(ct$tag, "_", reg$label)
    out_pdf <- file.path(OUT_DIR, paste0("final_lowres_", tag, ".pdf"))
    message("Plotting: ", tag)

    prof1 <- computeMethylationProfile(
      methylationData = meth_data[[ct$g1]],
      region          = reg$region,
      windowSize      = reg$window,
      context         = "CG")
    prof2 <- computeMethylationProfile(
      methylationData = meth_data[[ct$g2]],
      region          = reg$region,
      windowSize      = reg$window,
      context         = "CG")

    profiles <- GRangesList(
      setNames(list(prof1, prof2), c(ct$g1, ct$g2)))

    pdf(out_pdf, width = 10, height = 4)
    par(mar = c(4, 4, 3, 1) + 0.1)
    plotMethylationProfile(
      methylationProfiles = profiles,
      autoscale  = FALSE,
      title      = sprintf("CpG methylation %s (%d kb bins)\n%s vs %s",
                           reg$label, reg$window/1000, ct$g1, ct$g2),
      col        = unname(group_cols[c(ct$g1, ct$g2)]),
      pch        = c(15, 16),
      lty        = c(1, 2))
    dev.off()

    message("  Saved: ", out_pdf)
  }
}
message("All done.")
