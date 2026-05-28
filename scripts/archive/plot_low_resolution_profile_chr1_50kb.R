## Low-resolution methylation profile on chr1 (50 kb windows, one condition per plot)

## --- SETUP SECTION ---
## Copy the top part of your existing low-res script here:
## (library(DMRcaller), loading methylation objects, defining regions_chr1, etc.)
##
## For example, it might look like:
##
 library(DMRcaller)
## load("path/to/your_methylation_objects_chr1.RData")  # <- adjust to your file
## # objects created: methylationDataList, regions_chr1
##
## Make sure that after this section you have:
##   - methylationDataList[["Scramble_CTRL"]]
##   - methylationDataList[["ASO_VPA"]]
##   - regions_chr1 (GRanges covering chr1)

## --- 50 kb profiles ---
profile_ctrl <- computeMethylationProfile(
  methylationDataList[["Scramble_CTRL"]],
  regions_chr1,
  windowSize = 50000,
  context    = "CG"
)

profile_vpa <- computeMethylationProfile(
  methylationDataList[["ASO_VPA"]],
  regions_chr1,
  windowSize = 50000,
  context    = "CG"
)

profiles_ctrl <- GRangesList("Scramble_CTRL" = profile_ctrl)
profiles_vpa  <- GRangesList("ASO_VPA"       = profile_vpa)

dir.create("results/qc/dmrcaller", showWarnings = FALSE, recursive = TRUE)

## Scramble_CTRL only
pdf("results/qc/dmrcaller/low_resolution_profile_chr1_Scramble_CTRL_50kb.pdf",
    width = 10, height = 4)
par(mar = c(4, 4, 3, 1) + 0.1)
plotMethylationProfile(
  profiles_ctrl,
  autoscale = FALSE,
  title = "CG methylation on chr1 – Scramble_CTRL (50 kb windows)"
)
dev.off()

## ASO_VPA only
pdf("results/qc/dmrcaller/low_resolution_profile_chr1_ASO_VPA_50kb.pdf",
    width = 10, height = 4)
par(mar = c(4, 4, 3, 1) + 0.1)
plotMethylationProfile(
  profiles_vpa,
  autoscale = FALSE,
  title = "CG methylation on chr1 – ASO_VPA (50 kb windows)"
)
dev.off()
