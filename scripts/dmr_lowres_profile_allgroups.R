#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

OUT_DIR <- "results/qc/for_meeting"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
BS_DIR <- "results/alignments/bs"

group_cols <- c(
  ASO_VPA       = "#D94F3D",
  ASO_CTRL      = "#2E9B6F",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#1D6FA4"
)

message("Reading all groups...")
meth_data <- list(
  ASO_VPA       = readBismarkPool(Sys.glob(file.path(BS_DIR, "ASO_VPA_*_CX_report.txt.CpG_report.txt.gz"))),
  ASO_CTRL      = readBismarkPool(Sys.glob(file.path(BS_DIR, "ASO_CTRL_*_CX_report.txt.CpG_report.txt.gz"))),
  Scramble_VPA  = readBismarkPool(Sys.glob(file.path(BS_DIR, "Scramble_VPA_*_CX_report.txt.CpG_report.txt.gz"))),
  Scramble_CTRL = readBismarkPool(Sys.glob(file.path(BS_DIR, "Scramble_CTRL_*_CX_report.txt.CpG_report.txt.gz")))
)
message("All groups loaded.")

regions_list <- list(
  list(region = GRanges("chr1",  IRanges(1, 2.5e8)),
       window = 1000000, label = "chr1_1Mb"),
  list(region = GRanges("chr13", IRanges(60e6, 80e6)),
       window = 50000,   label = "chr13_60_80Mb_50kb")
)

contrasts <- list(
  list(g1="ASO_VPA",     g2="Scramble_CTRL", tag="combined_effect"),
  list(g1="Scramble_VPA",g2="Scramble_CTRL", tag="VPA_effect"),
  list(g1="ASO_VPA",     g2="ASO_CTRL",      tag="ASO_effect"),
  list(g1="ASO_CTRL",    g2="Scramble_CTRL", tag="baseline")
)

for (reg in regions_list) {
  for (ct in contrasts) {
    out_pdf <- file.path(OUT_DIR,
      paste0("final_lowres_", ct$tag, "_", reg$label, ".pdf"))
    message("Plotting: ", ct$tag, " — ", reg$label)

    prof_list <- lapply(c(ct$g1, ct$g2), function(g) {
      prof <- computeMethylationProfile(
        methylationData = meth_data[[g]],
        region          = reg$region,
        windowSize      = reg$window,
        context         = "CG")
      df <- as.data.frame(prof)
      df$meth <- df$sumReadsM / df$sumReadsN
      df$pos  <- (df$start + df$end) / 2
      df
    })
    names(prof_list) <- c(ct$g1, ct$g2)

    x_range <- range(do.call(c, lapply(prof_list, function(d) d$pos)), na.rm=TRUE)

    pdf(out_pdf, width = 10, height = 4)
    par(mar = c(4, 4, 3, 1) + 0.1)
    plot(NULL, xlim = x_range, ylim = c(0, 1),
         xlab = "genomic coordinate", ylab = "methylation",
         main = sprintf("CpG methylation %s (%d kb bins)\n%s vs %s",
                        reg$label, reg$window/1000, ct$g1, ct$g2))
    for (g in c(ct$g1, ct$g2)) {
      d <- prof_list[[g]]
      lines(d$pos, d$meth, col = group_cols[g], lwd = 1.5)
    }
    legend("topright", legend = c(ct$g1, ct$g2),
           col = group_cols[c(ct$g1, ct$g2)],
           lwd = 1.5, bty = "n", cex = 0.85)
    dev.off()
    message("  Saved: ", out_pdf)
  }
}
message("All done.")
