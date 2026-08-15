.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

COV_DIR <- "results/alignments/bs/by_chr"

all_groups <- list(
  ASO_VPA       = Sys.glob(file.path(COV_DIR, "ASO_VPA_*_chr1.CpG_report.txt.gz")),
  ASO_CTRL      = Sys.glob(file.path(COV_DIR, "ASO_CTRL_*_chr1.CpG_report.txt.gz")),
  Scramble_VPA  = Sys.glob(file.path(COV_DIR, "Scramble_VPA_*_chr1.CpG_report.txt.gz")),
  Scramble_CTRL = Sys.glob(file.path(COV_DIR, "Scramble_CTRL_*_chr1.CpG_report.txt.gz"))
)

contrasts <- list(
  list(g1="ASO_VPA",      g2="ASO_CTRL",      tag="ASO_effect"),
  list(g1="Scramble_VPA", g2="Scramble_CTRL",  tag="VPA_effect"),
  list(g1="ASO_VPA",      g2="Scramble_CTRL",  tag="combined_effect")
)

region      <- GRanges(seqnames=Rle("chr1"), ranges=IRanges(1, 10e6))
window_size <- 50000
out_dir     <- "results/qc/dmrcaller"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Load all groups once
message("Loading all groups...")
profiles_cache <- list()
for (nm in names(all_groups)) {
  message("  Reading ", nm, "...")
  meth <- readBismarkPool(all_groups[[nm]])
  prof <- computeMethylationProfile(
    methylationData=meth, region=region,
    windowSize=window_size, context="CG")
  profiles_cache[[nm]] <- prof
  rm(meth); gc()
}

colours <- c("#E69F00", "#56B4E9")

for (ct in contrasts) {
  g1 <- ct$g1; g2 <- ct$g2; tag <- ct$tag
  message("Plotting: ", g1, " vs ", g2)

  profiles_gr <- GRangesList(setNames(
    list(profiles_cache[[g1]], profiles_cache[[g2]]),
    c(g1, g2)))

  out_pdf <- file.path(out_dir, paste0("low_resolution_", tag, "_chr1_50kb.pdf"))
  pdf(out_pdf, width=10, height=5)
  par(mar=c(4,4,3,1)+0.1)
  plotMethylationProfile(
    methylationProfiles = profiles_gr,
    autoscale = FALSE,
    labels    = NULL,
    title     = paste0("CG methylation chr1:1-10Mb (50kb bins)\n", g1, " vs ", g2),
    col       = colours,
    pch       = c(15, 16),
    lty       = c(1, 2))
  dev.off()
  message("Saved: ", out_pdf)
}
message("Done.")
