.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

OUT_DIR    <- "results/lowres_profiles"
METH_CACHE <- "results/dmr/meth_pooled_cache.rds"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

GROUP_COLS <- c(
  ASO_VPA       = "#B2182B",
  ASO_CTRL      = "#1B4F8A",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#6B7280"
)
GROUP_LTY <- c(
  ASO_VPA       = 1,
  ASO_CTRL      = 1,
  Scramble_VPA  = 2,
  Scramble_CTRL = 2
)
GROUPS <- c("ASO_VPA", "ASO_CTRL", "Scramble_VPA", "Scramble_CTRL")

# TIER 1: chromosome overviews (250 kb bins) -- chr13 swapped for chr17 CACNG cluster
# chrX=3.97/Mb (top), chr8=1.74/Mb (CHRNB3), chr17=CACNG cluster (calcium/AMPA signature),
# chr14=1.02/Mb (MTA1-DT), chr5=0.84/Mb (SMN2, primary biological target)
REGIONS <- list(
  list(region=GRanges("chrX",  IRanges(1, 156040895)),
       window=250000, label="chrX_250kb",
       title="CpG methylation chrX (250 kb bins): top DMR density 3.97/Mb"),
  list(region=GRanges("chr8",  IRanges(1, 145138636)),
       window=250000, label="chr8_250kb",
       title="CpG methylation chr8 (250 kb bins): 1.74 DMRs/Mb, CHRNB3 locus"),
  list(region=GRanges("chr17", IRanges(1, 83257441)),
       window=250000, label="chr17_250kb",
       title="CpG methylation chr17 (250 kb bins): CACNG/PRKCA calcium signalling cluster"),
  list(region=GRanges("chr14", IRanges(1, 107043718)),
       window=250000, label="chr14_250kb",
       title="CpG methylation chr14 (250 kb bins): 1.02 DMRs/Mb, MTA1-DT locus"),
  list(region=GRanges("chr5",  IRanges(1, 181538259)),
       window=250000, label="chr5_250kb",
       title="CpG methylation chr5 (250 kb bins): SMN2 primary target"),

  # TIER 2: regional zooms (25 kb bins) -- centred on known hotspots
  list(region=GRanges("chrX",  IRanges(5000000, 30000000)),
       window=25000,  label="chrX_5_30Mb_25kb",
       title="CpG methylation chrX:5-30 Mb (25 kb bins): ASO DMR cluster"),
  list(region=GRanges("chr8",  IRanges(38000000, 55000000)),
       window=25000,  label="chr8_CHRNB3_25kb",
       title="CpG methylation chr8:38-55 Mb (25 kb bins): CHRNB3 region"),
  list(region=GRanges("chr17", IRanges(66000000, 67500000)),
       window=25000,  label="chr17_CACNG_25kb",
       title="CpG methylation chr17:66-67.5 Mb (25 kb bins): PRKCA-CACNG5-CACNG4-CACNG1"),
  list(region=GRanges("chr14", IRanges(100000000, 107000000)),
       window=25000,  label="chr14_MTA1DT_25kb",
       title="CpG methylation chr14:100-107 Mb (25 kb bins): MTA1-DT region"),
  list(region=GRanges("chr5",  IRanges(67000000, 73000000)),
       window=25000,  label="chr5_SMN2_25kb",
       title="CpG methylation chr5:67-73 Mb (25 kb bins): SMN2 region")
)

message("loading pooled methylation cache...")
if (!file.exists(METH_CACHE)) stop("cache not found: ", METH_CACHE)
meth_pooled <- readRDS(METH_CACHE)
message("  loaded: ", paste(names(meth_pooled), collapse=", "))

for (reg in REGIONS) {
  message("\nregion: ", reg$label)

  prof_list <- list()
  for (g in GROUPS) {
    message("  computing: ", g)
    prof        <- computeMethylationProfile(meth_pooled[[g]], reg$region,
                                             windowSize=reg$window, context="CG")
    df          <- as.data.frame(prof)
    df$meth     <- df$sumReadsM / df$sumReadsN
    df$pos      <- (df$start + df$end) / 2
    df$group    <- g
    prof_list[[g]] <- df[!is.na(df$meth) & df$sumReadsN >= 3, ]
  }
  x_range <- range(unlist(lapply(prof_list, function(d) d$pos)), na.rm=TRUE)

  out_pdf <- file.path(OUT_DIR, paste0("lowres_allgroups_", reg$label, ".pdf"))
  pdf(out_pdf, width=10, height=4)
  par(mar=c(4,4,3,9)+0.1, bg="white", col.axis="black",
      col.lab="black", col.main="black", fg="black")
  plot(NULL, xlim=x_range, ylim=c(0,1), xaxt="n",
       xlab="genomic coordinate (Mb)", ylab="CpG methylation proportion",
       main=reg$title, cex.main=0.9, font.main=2)
  axis_at <- pretty(x_range, n=8)
  axis(1, at=axis_at, labels=sprintf("%g", axis_at/1e6))
  for (g in GROUPS) {
    d <- prof_list[[g]]
    if (nrow(d) > 0) lines(d$pos, d$meth, col=GROUP_COLS[g], lwd=2.0, lty=GROUP_LTY[g])
  }
  legend("topright", legend=GROUPS, col=GROUP_COLS[GROUPS], lty=GROUP_LTY[GROUPS],
         lwd=1.8, bty="n", cex=0.85, xpd=TRUE, inset=c(-0.22,0))
  dev.off()
  message("  saved: ", basename(out_pdf))
}

message("\nDone.")
