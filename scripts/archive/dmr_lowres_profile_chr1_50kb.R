#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

group_files <- list(
  ASO_VPA = Sys.glob("results/alignments/bs/ASO_VPA_*_CX_report.txt.CpG_report.txt.gz"),
  ASO_CTRL = Sys.glob("results/alignments/bs/ASO_CTRL_*_CX_report.txt.CpG_report.txt.gz"),
  Scramble_VPA = Sys.glob("results/alignments/bs/Scramble_VPA_*_CX_report.txt.CpG_report.txt.gz"),
  Scramble_CTRL = Sys.glob("results/alignments/bs/Scramble_CTRL_*_CX_report.txt.CpG_report.txt.gz")
)

regions     <- GRanges(seqnames=Rle("chr1"), ranges=IRanges(1, 3e8))
window_size <- 50000
out_dir     <- "results/qc/dmrcaller"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

for (nm in names(group_files)) {
  out_pdf <- sprintf("%s/low_resolution_profile_chr1_%s_%dkb.pdf",
                     out_dir, nm, window_size/1000)

  # Skip if already done
  if (file.exists(out_pdf)) {
    message("Skipping ", nm, " — already exists")
    next
  }

  fs <- group_files[[nm]]
  if (length(fs) == 0) { message("Skipping ", nm, " (no files)"); next }

  message("Reading & pooling group: ", nm)
  meth <- readBismarkPool(fs)
  message("Computing profile for ", nm)
  prof <- computeMethylationProfile(
    methylationData = meth,
    region          = regions,
    windowSize      = window_size,
    context         = "CG")

  profiles_one <- GRangesList(setNames(list(prof), nm))

  pdf(out_pdf, width=8, height=4)
  par(mar=c(4,4,3,1)+0.1)
  plotMethylationProfile(
    methylationProfiles = profiles_one,
    autoscale = FALSE,
    labels    = NULL,
    title     = sprintf("CG methylation on chr1 - %s (%d kb bins)", nm, window_size/1000),
    col       = "black", pch=16, lty=1)
  dev.off()
  message("Saved: ", out_pdf)
  rm(meth, prof, profiles_one); gc()
}
message("Done.")
