#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

# Low-resolution methylation profiles across all 4 groups.
# Uses the pooled methylation cache from 07b_dmr_plots.R to avoid re-reading raw files.
# Generates: one all-groups plot per region + one pairwise contrast plot per contrast per region.
# Regions: chr1 overview, chrX overview + zoom (ASO DMR hotspot), chr5 + SMN2 zoom.

METH_CACHE <- "results/dmr/meth_pooled_cache.rds"
OUT_DIR    <- "results/lowres_profiles"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

GROUP_COLS <- c(
  ASO_VPA       = "#B2182B",
  ASO_CTRL      = "#1B4F8A",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#6B7280"
)
GROUPS <- c("ASO_VPA", "ASO_CTRL", "Scramble_VPA", "Scramble_CTRL")

REGIONS <- list(
  list(region=GRanges("chr1",  IRanges(1, 248956422)),
       window=1000000, label="chr1_1Mb",
       title="CpG methylation chr1 (1 Mb bins)"),
  list(region=GRanges("chrX",  IRanges(1, 155270560)),
       window=500000,  label="chrX_500kb",
       title="CpG methylation chrX (500 kb bins) — ASO hotspot (620 DMRs)"),
  list(region=GRanges("chrX",  IRanges(10000000, 20000000)),
       window=50000,   label="chrX_10_20Mb_50kb",
       title="CpG methylation chrX:10-20 Mb (50 kb bins) — highest ASO density"),
  list(region=GRanges("chr5",  IRanges(1, 181538259)),
       window=500000,  label="chr5_500kb",
       title="CpG methylation chr5 (500 kb bins)"),
  list(region=GRanges("chr5",  IRanges(69000000, 71500000)),
       window=10000,   label="chr5_SMN2_10kb",
       title="CpG methylation chr5:69-71.5 Mb (10 kb bins) — SMN2 locus")
)

CONTRASTS <- list(
  list(g1="ASO_CTRL",     g2="Scramble_CTRL", tag="ASO_effect",
       title="ASO effect (nusinersen vs scramble CTRL)"),
  list(g1="Scramble_VPA", g2="Scramble_CTRL", tag="VPA_effect",
       title="VPA effect (HDAC inhibitor vs CTRL)"),
  list(g1="ASO_VPA",      g2="Scramble_CTRL", tag="combination",
       title="Combination effect (ASO+VPA vs CTRL)"),
  list(g1="ASO_VPA",      g2="ASO_CTRL",      tag="VPA_on_ASO",
       title="VPA effect on ASO background"),
  list(g1="ASO_VPA",      g2="Scramble_VPA",  tag="ASO_on_VPA",
       title="ASO effect on VPA background")
)

message("loading pooled methylation cache...")
if (!file.exists(METH_CACHE)) stop("cache not found: ", METH_CACHE)
meth_pooled <- readRDS(METH_CACHE)
message("  loaded: ", paste(names(meth_pooled), collapse=", "))

for (reg in REGIONS) {
  message("\nregion: ", reg$label)

  # compute profile for each group
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

  # all-groups plot
  out_pdf <- file.path(OUT_DIR, paste0("lowres_allgroups_", reg$label, ".pdf"))
  pdf(out_pdf, width=10, height=4)
  par(mar=c(4,4,3,8)+0.1, bg="white", col.axis="black",
      col.lab="black", col.main="black", fg="black")
  plot(NULL, xlim=x_range, ylim=c(0,1),
       xlab="genomic coordinate", ylab="CpG methylation proportion",
       main=reg$title, cex.main=0.9, font.main=2)
  for (g in GROUPS) {
    d <- prof_list[[g]]
    if (nrow(d) > 0) lines(d$pos, d$meth, col=GROUP_COLS[g], lwd=1.5)
  }
  legend("topright", legend=GROUPS, col=GROUP_COLS[GROUPS],
         lwd=1.5, bty="n", cex=0.85, xpd=TRUE, inset=c(-0.15,0))
  dev.off()
  message("  saved: ", basename(out_pdf))

  # pairwise contrast plots
  for (ct in CONTRASTS) {
    out_pdf <- file.path(OUT_DIR,
      paste0("lowres_", ct$tag, "_", reg$label, ".pdf"))
    pdf(out_pdf, width=10, height=4)
    par(mar=c(4,4,3,2)+0.1, bg="white", col.axis="black",
        col.lab="black", col.main="black", fg="black")
    plot(NULL, xlim=x_range, ylim=c(0,1),
         xlab="genomic coordinate", ylab="CpG methylation proportion",
         main=sprintf("%s\n%s (%d kb bins)", ct$title, reg$label, reg$window/1000),
         cex.main=0.9, font.main=2)
    for (g in c(ct$g1, ct$g2)) {
      d <- prof_list[[g]]
      if (nrow(d) > 0) lines(d$pos, d$meth, col=GROUP_COLS[g], lwd=1.5)
    }
    legend("topright", legend=c(ct$g1, ct$g2),
           col=GROUP_COLS[c(ct$g1, ct$g2)],
           lwd=1.5, bty="n", cex=0.85)
    dev.off()
  }
}
message("\ndone. outputs in: ", OUT_DIR)
