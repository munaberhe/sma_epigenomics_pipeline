#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

# 1) Group-specific CX report files (full CX reports)[web:212]
group_files <- list(
  ASO_VPA       = Sys.glob("results/alignments/bs/ASO_VPA_*bismark.deduplicated.CX_report.txt*"),
  ASO_CTRL      = Sys.glob("results/alignments/bs/ASO_CTRL_*bismark.deduplicated.CX_report.txt*"),
  Scramble_VPA  = Sys.glob("results/alignments/bs/Scramble_VPA_*bismark.deduplicated.CX_report.txt*"),
  Scramble_CTRL = Sys.glob("results/alignments/bs/Scramble_CTRL_*bismark.deduplicated.CX_report.txt*")
)

print(group_files)

meth_list <- lapply(group_files, function(fs) {
  if (length(fs) == 0) stop("No CX reports found for one group")
  message("Reading & pooling: ", paste(fs, collapse = ", "))
  # Reads (gzipped) Bismark CX reports into a methylationData object[web:212]
  readBismarkPool(fs)
})

# 2) Define chr1 region (similar to DMRcaller vignette, but for chr1)[web:214]
regions <- GRanges(seqnames = Rle("chr1"),
                   ranges   = IRanges(1, 1e8))  # large upper bound for chr1

# 3) Compute low-res profiles in 10 kb windows, CG context[web:211]
window_size <- 10000
profilesCG <- GRangesList(lapply(names(meth_list), function(nm) {
  message("Computing profile for ", nm)
  computeMethylationProfile(
    methylationData = meth_list[[nm]],
    region          = regions,
    windowSize      = window_size,
    context         = "CG"
  )
}))
names(profilesCG) <- names(meth_list)

# 4) Plot low-resolution profile like vignette's example[web:216]
out_pdf <- "results/qc/dmrcaller/low_resolution_profile_chr1_DMRcaller.pdf"
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)

pdf(out_pdf, width = 8, height = 4)
par(mar = c(4, 4, 3, 1) + 0.1)
par(mfrow = c(1, 1))

plotMethylationProfile(
  methylationProfiles = profilesCG,
  autoscale = FALSE,
  labels    = NULL,
  title     = paste0("CG methylation on Chromosome 1 (", window_size/1000, " kb bins)"),
  col       = c("#D55E00","#E69F00","#0072B2","#56B4E9"),
  pch       = c(1, 0, 16, 2),
  lty       = c(4, 1, 3, 2)
)
dev.off()
message("Saved: ", out_pdf)
