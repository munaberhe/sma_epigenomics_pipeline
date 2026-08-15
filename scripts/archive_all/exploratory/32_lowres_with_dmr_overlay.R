.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

OUT_DIR    <- "results/lowres_profiles"
METH_CACHE <- "results/dmr/meth_pooled_cache.rds"
DMR_DIR    <- "results/dmr_annotation"
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

# Regional zooms (25 kb bins) with real DMR overlay.
# Each region pulls significant DMRs from the same annotated CSV used
# throughout this analysis -- combination contrast (ASO_VPA_vs_Scramble_CTRL)
# shown as red ticks, ASO-alone contrast (ASO_CTRL_vs_Scramble_CTRL) as navy ticks.
REGIONS <- list(
  list(region=GRanges("chrX",  IRanges(5000000, 30000000)),
       window=25000,  label="chrX_5_30Mb_25kb_dmr",
       title="CpG methylation chrX:5-30 Mb (25 kb bins): ASO DMR cluster"),
  list(region=GRanges("chr8",  IRanges(38000000, 55000000)),
       window=25000,  label="chr8_CHRNB3_25kb_dmr",
       title="CpG methylation chr8:38-55 Mb (25 kb bins): CHRNB3 region"),
  list(region=GRanges("chr17", IRanges(66000000, 67500000)),
       window=25000,  label="chr17_CACNG_25kb_dmr",
       title="CpG methylation chr17:66-67.5 Mb (25 kb bins): PRKCA-CACNG5-CACNG4-CACNG1"),
  list(region=GRanges("chr14", IRanges(100000000, 107000000)),
       window=25000,  label="chr14_MTA1DT_25kb_dmr",
       title="CpG methylation chr14:100-107 Mb (25 kb bins): MTA1-DT region"),
  list(region=GRanges("chr5",  IRanges(67000000, 73000000)),
       window=25000,  label="chr5_SMN2_25kb_dmr",
       title="CpG methylation chr5:67-73 Mb (25 kb bins): SMN2 region")
)

# load DMR positions for a region from the annotated CSV
load_dmr_positions <- function(csv_path, region) {
  if (!file.exists(csv_path)) {
    message("    DMR CSV not found: ", csv_path)
    return(numeric(0))
  }
  df <- read.csv(csv_path)
  chr <- as.character(seqnames(region))
  rstart <- start(region); rend <- end(region)
  sub <- df[df$seqnames == chr & df$end >= rstart & df$start <= rend, ]
  (sub$start + sub$end) / 2
}

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

  # load real DMR positions for this region
  combo_dmrs <- load_dmr_positions(
    file.path(DMR_DIR, "ASO_VPA_vs_Scramble_CTRL_annotated.csv"), reg$region)
  aso_dmrs <- load_dmr_positions(
    file.path(DMR_DIR, "ASO_CTRL_vs_Scramble_CTRL_annotated.csv"), reg$region)
  message("  combination DMRs in region: ", length(combo_dmrs),
          "  ASO-alone DMRs in region: ", length(aso_dmrs))

  out_pdf <- file.path(OUT_DIR, paste0("lowres_allgroups_", reg$label, ".pdf"))
  pdf(out_pdf, width=10, height=5)
  par(mar=c(4,4,5,9)+0.1, bg="white", col.axis="black",
      col.lab="black", col.main="black", fg="black")
  plot(NULL, xlim=x_range, ylim=c(-0.08,1), xaxt="n",
       xlab="genomic coordinate (Mb)", ylab="CpG methylation proportion",
       main=reg$title, cex.main=0.9, font.main=2)
  axis_at <- pretty(x_range, n=8)
  axis(1, at=axis_at, labels=sprintf("%g", axis_at/1e6))

  # DMR tick marks below the x-axis line, two rows
  if (length(combo_dmrs) > 0)
    points(combo_dmrs, rep(-0.03, length(combo_dmrs)), pch="|",
           col=GROUP_COLS["ASO_VPA"], cex=0.9)
  if (length(aso_dmrs) > 0)
    points(aso_dmrs, rep(-0.06, length(aso_dmrs)), pch="|",
           col=GROUP_COLS["ASO_CTRL"], cex=0.9)
  mtext("DMRs:", side=2, at=-0.045, las=1, cex=0.6, line=2.2)

  for (g in GROUPS) {
    d <- prof_list[[g]]
    if (nrow(d) > 0) lines(d$pos, d$meth, col=GROUP_COLS[g], lwd=2.0, lty=GROUP_LTY[g])
  }
  legend("topright", legend=c(GROUPS, "Combo DMR", "ASO-alone DMR"),
         col=c(GROUP_COLS[GROUPS], GROUP_COLS["ASO_VPA"], GROUP_COLS["ASO_CTRL"]),
         lty=c(GROUP_LTY[GROUPS], NA, NA),
         pch=c(NA,NA,NA,NA,"|","|"),
         lwd=1.8, bty="n", cex=0.8, xpd=TRUE, inset=c(-0.27,0))
  dev.off()
  message("  saved: ", basename(out_pdf))
}

message("\nDone.")
