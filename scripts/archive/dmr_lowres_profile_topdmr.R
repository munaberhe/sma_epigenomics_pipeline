.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
})

COV_DIR <- "results/alignments/bs/by_chr"
out_dir <- "results/qc/dmrcaller"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

all_groups <- list(
  ASO_VPA       = Sys.glob(file.path(COV_DIR, "ASO_VPA_*_chr1.CpG_report.txt.gz")),
  ASO_CTRL      = Sys.glob(file.path(COV_DIR, "ASO_CTRL_*_chr1.CpG_report.txt.gz")),
  Scramble_VPA  = Sys.glob(file.path(COV_DIR, "Scramble_VPA_*_chr1.CpG_report.txt.gz")),
  Scramble_CTRL = Sys.glob(file.path(COV_DIR, "Scramble_CTRL_*_chr1.CpG_report.txt.gz"))
)

contrasts <- list(
  list(g1="ASO_VPA",      g2="ASO_CTRL",      tag="topdmr_chr1_50kb",
       colours=c("ASO_VPA"="#D55E00","ASO_CTRL"="#0072B2"),
       title="CG methylation chr1:224-236Mb (50kb bins)\nASO_VPA vs ASO_CTRL"),
  list(g1="Scramble_VPA", g2="Scramble_CTRL",  tag="topdmr_VPA_effect_chr1_50kb",
       colours=c("Scramble_VPA"="#CC79A7","Scramble_CTRL"="#009E73"),
       title="CG methylation chr1:224-236Mb (50kb bins)\nScramble_VPA vs Scramble_CTRL"),
  list(g1="ASO_VPA",      g2="Scramble_CTRL",  tag="topdmr_combined_effect_chr1_50kb",
       colours=c("ASO_VPA"="#E69F00","Scramble_CTRL"="#56B4E9"),
       title="CG methylation chr1:224-236Mb (50kb bins)\nASO_VPA vs Scramble_CTRL")
)

region      <- GRanges(seqnames=Rle("chr1"), ranges=IRanges(224e6, 236e6))
window_size <- 50000

message("Loading all groups...")
profiles_cache <- list()
for (nm in names(all_groups)) {
  message("  Reading ", nm, "...")
  meth <- readBismarkPool(all_groups[[nm]])
  prof <- computeMethylationProfile(
    methylationData=meth, region=region,
    windowSize=window_size, context="CG")
  d <- as.data.frame(prof)
  d$pos_mb <- (d$start + d$end) / 2 / 1e6
  d$methylation <- d$sumReadsM / d$sumReadsN
  d$group <- nm
  profiles_cache[[nm]] <- d
  rm(meth, prof); gc()
}

for (ct in contrasts) {
  g1 <- ct$g1; g2 <- ct$g2
  message("Plotting: ", g1, " vs ", g2)

  dat <- rbind(profiles_cache[[g1]], profiles_cache[[g2]])
  dat$group <- factor(dat$group, levels=c(g1, g2))

  p <- ggplot(dat, aes(x=pos_mb, y=methylation, colour=group, group=group)) +
    geom_line(linewidth=0.9) +
    geom_point(size=1.5) +
    scale_colour_manual(values=ct$colours, name=NULL) +
    scale_x_continuous(breaks=seq(224, 236, 2)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
    labs(title=ct$title, x="chr1 position (Mb)", y="Mean CpG methylation") +
    theme_classic(base_size=12) +
    theme(legend.position="bottom",
          legend.title=element_blank(),
          plot.title=element_text(face="bold", size=11))

  out_pdf <- file.path(out_dir, paste0("low_resolution_", ct$tag, ".pdf"))
  ggsave(out_pdf, p, width=10, height=5)
  message("Saved: ", out_pdf)
}
message("Done.")
