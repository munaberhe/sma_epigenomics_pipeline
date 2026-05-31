#!/usr/bin/env Rscript
# 15_lowres_methylation_profile.R
# Low-resolution methylation profiles in the style approved by Radu (May 2026)
# Uses computeMethylationProfile() from DMRcaller on the pooled cache
# so we don't need to re-read raw files -- saves ~25 minutes.
# Plots all 4 groups as coloured lines on the same axes.
# Regions: chr1 (1Mb bins, genome overview), chrX (500kb + 50kb zoom), chr5 (500kb + SMN2 10kb zoom)
# chrX chosen because ASO_CTRL has 620 DMRs there = 18% of all ASO DMRs,
# disproportionate to chrX size (~5% of genome) -- genuine hotspot.

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
.libPaths(c("~/R/library", .libPaths()))

METH_CACHE <- "results/dmr/meth_pooled_cache.rds"
OUT_DIR    <- "results/lowres_profiles"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# Consistent colour scheme across all thesis figures
GROUP_COLS <- c(
  ASO_VPA       = "#D94F3D",
  ASO_CTRL      = "#2E9B6F",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#1D6FA4"
)

# Regions to plot -- chr1 for genome-wide overview, chrX for ASO hotspot.
# chrX 10-20Mb chosen because it has the highest ASO DMR density (93 DMRs)
REGIONS <- list(
  list(
    region = GRanges("chr1", IRanges(1, 248956422)),
    window = 1000000,
    label  = "chr1_1Mb",
    title  = "CpG methylation chr1 (1 Mb bins)\nAll groups"
  ),
  list(
    region = GRanges("chrX", IRanges(1, 155270560)),
    window = 500000,
    label  = "chrX_500kb",
    title  = "CpG methylation chrX (500 kb bins)\nAll groups -- ASO DMR hotspot (620 DMRs, 18% of total)"
  ),
  list(
    region = GRanges("chrX", IRanges(10000000, 20000000)),
    window = 50000,
    label  = "chrX_10_20Mb_50kb",
    title  = "CpG methylation chrX:10-20 Mb (50 kb bins)\nHighest ASO DMR density window (93 DMRs)"
    ),
  list(
    region = GRanges("chr5", IRanges(1, 181538259)),
    window = 500000,
    label  = "chr5_500kb",
    title  = "CpG methylation chr5 (500 kb bins)\nSMN2 locus chromosome overview"
  ),
  list(
    region = GRanges("chr5", IRanges(69000000, 71500000)),
    window = 10000,
    label  = "chr5_SMN2_10kb",
    title  = "CpG methylation chr5:69-71.5 Mb (10 kb bins)\nSMN2 locus zoom -- all 4 conditions"
  )
)

# Load pooled methylation cache
message("Loading pooled methylation cache...")
if (!file.exists(METH_CACHE)) stop("Cache not found: ", METH_CACHE)
meth_pooled <- readRDS(METH_CACHE)
message("  Loaded: ", paste(names(meth_pooled), collapse=", "))

# All four groups plotted together -- same order as the Radu-approved slide
GROUPS <- c("ASO_VPA", "ASO_CTRL", "Scramble_VPA", "Scramble_CTRL")

# Generate one PDF per region
for (reg in REGIONS) {
  message("\nPlotting: ", reg$label)

  out_pdf <- file.path(OUT_DIR, paste0("lowres_allgroups_", reg$label, ".pdf"))

  # Compute methylation profile for each group using DMRcaller
  # computeMethylationProfile() tiles the region into bins and calculates
  # the mean CpG methylation proportion in each bin
  prof_list <- list()
  for (g in GROUPS) {
    message("  Computing profile: ", g)
    prof <- computeMethylationProfile(
      methylationData = meth_pooled[[g]],
      region          = reg$region,
      windowSize      = reg$window,
      context         = "CG"
    )
    df       <- as.data.frame(prof)
    df$meth  <- df$sumReadsM / df$sumReadsN
    df$pos   <- (df$start + df$end) / 2
    df$group <- g
    # Only keep bins with at least 3 reads to avoid noise from low-coverage bins
    prof_list[[g]] <- df[!is.na(df$meth) & df$sumReadsN >= 3, ]
  }

  # Get x range from data
  x_range <- range(unlist(lapply(prof_list, function(d) d$pos)), na.rm=TRUE)

  pdf(out_pdf, width=10, height=4)
  par(mar=c(4, 4, 3, 8) + 0.1, bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")

  # Empty plot to set axes
  plot(NULL, xlim=x_range, ylim=c(0, 1),
       xlab="genomic coordinate",
       ylab="CpG methylation proportion",
       main=reg$title,
       cex.main=0.9, font.main=2)

  # Draw one line per group
  for (g in GROUPS) {
    d <- prof_list[[g]]
    if (nrow(d) == 0) next
    lines(d$pos, d$meth, col=GROUP_COLS[g], lwd=1.5)
  }

  # Legend outside plot area on the right
  legend("topright", legend=GROUPS,
         col=GROUP_COLS[GROUPS],
         lwd=1.5, bty="n", cex=0.85,
         xpd=TRUE, inset=c(-0.15, 0))

  dev.off()
  message("  Saved: ", out_pdf)
}

# Also generate pairwise contrast plots for each region
# These show each treatment contrast directly -- cleaner for thesis figures
CONTRASTS <- list(
  list(g1="ASO_CTRL",    g2="Scramble_CTRL", tag="ASO_effect",      title="ASO effect (nusinersen vs scramble CTRL)"),
  list(g1="Scramble_VPA",g2="Scramble_CTRL", tag="VPA_effect",      title="VPA effect (HDAC inhibitor vs CTRL)"),
  list(g1="ASO_VPA",     g2="Scramble_CTRL", tag="combination",     title="Combination effect (ASO+VPA vs CTRL)"),
  list(g1="ASO_VPA",     g2="ASO_CTRL",      tag="VPA_on_ASO",      title="VPA effect on ASO background (ASO+VPA vs ASO alone)"),
  list(g1="ASO_VPA",     g2="Scramble_VPA",  tag="ASO_on_VPA",      title="ASO effect on VPA background (ASO+VPA vs VPA alone)")
)

for (reg in REGIONS) {
  for (ct in CONTRASTS) {
    message("Plotting: ", ct$tag, " -- ", reg$label)

    out_pdf <- file.path(OUT_DIR,
      paste0("lowres_", ct$tag, "_", reg$label, ".pdf"))

    x_range <- range(unlist(lapply(
      list(prof_list[[ct$g1]], prof_list[[ct$g2]]),
      function(d) d$pos)), na.rm=TRUE)

    pdf(out_pdf, width=10, height=4)
    par(mar=c(4, 4, 3, 2) + 0.1, bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")

    plot(NULL, xlim=x_range, ylim=c(0, 1),
         xlab="genomic coordinate",
         ylab="CpG methylation proportion",
         main=sprintf("%s\n%s (%d kb bins)",
                      ct$title, reg$label, reg$window/1000),
         cex.main=0.9, font.main=2)

    for (g in c(ct$g1, ct$g2)) {
      d <- prof_list[[g]]
      if (nrow(d) == 0) next
      lines(d$pos, d$meth, col=GROUP_COLS[g], lwd=1.5)
    }

    legend("topright",
           legend=c(ct$g1, ct$g2),
           col=GROUP_COLS[c(ct$g1, ct$g2)],
           lwd=1.5, bty="n", cex=0.85)

    dev.off()
    message("  Saved: ", out_pdf)
  }
}

message("\nAll done. Outputs in: ", OUT_DIR)
list.files(OUT_DIR, pattern="\\.pdf$")
