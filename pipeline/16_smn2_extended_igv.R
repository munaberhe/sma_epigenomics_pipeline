#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(karyoploteR)
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/figures/smn2_extended_igv"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)


# SMN2: chr5:70,049,638-70,078,522
SMN2_START <- 70049638
SMN2_END   <- 70078522
FLANK      <- 200000
REGION     <- GRanges("chr5", IRanges(SMN2_START - FLANK, SMN2_END + FLANK))


COND_COLS <- c(
  ASO_CTRL      = "#1F3A5F",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#C0392B",
  Scramble_VPA  = "#F0A500"
)
ENH_H9_COL    <- "#E69F00"
ENH_CCRE_COL  <- "#56B4E9"
CGI_COL       <- "#A8D5E2"
H3K27_CTRL    <- "#4E9EC7"
H3K27_VPA     <- "#8E44AD"
DMR_COL       <- "#E31A1C"


message("Loading annotation tracks...")

# CpG islands
cgi_raw <- read.table("data/reference/cpg_islands_hg38.bed", header=FALSE, sep="\t", stringsAsFactors=FALSE)
cgi <- GRanges(seqnames=cgi_raw[,1], ranges=IRanges(cgi_raw[,2], cgi_raw[,3]))
cgi <- subsetByOverlaps(cgi, REGION)
message("  CpG islands: ", length(cgi))

# H9 enhancers
h9_raw <- read.table(gzfile("data/reference/H9_predicted_non_promoter_non_fragments.bed.gz"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
h9_enh <- GRanges(seqnames=h9_raw$seqnames, ranges=IRanges(h9_raw$start, h9_raw$end))
h9_enh <- subsetByOverlaps(h9_enh, REGION)
message("  H9 enhancers: ", length(h9_enh))

# ENCODE cCREs
ccre_raw <- read.table("data/reference/encode_cCREs_hg38.bed", header=FALSE, sep="\t", stringsAsFactors=FALSE)
ccre <- GRanges(seqnames=ccre_raw[,1], ranges=IRanges(ccre_raw[,2], ccre_raw[,3]))
ccre <- subsetByOverlaps(ccre, REGION)
message("  ENCODE cCREs: ", length(ccre))

# H3K27ac peaks
h3k27_ctrl <- tryCatch({
  r1 <- import("data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep1.narrowPeak.gz",
               format="narrowPeak")
  r2 <- import("data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep2.narrowPeak.gz",
               format="narrowPeak")
  subsetByOverlaps(c(r1, r2), REGION)
}, error=function(e) { message("  H3K27ac CTRL not loaded: ", e$message); GRanges() })

h3k27_vpa <- tryCatch({
  r1 <- import("data/external/h3k27ac_gse246399/H3K27ac_VPA_Rep1.narrowPeak.gz",
               format="narrowPeak")
  subsetByOverlaps(r1, REGION)
}, error=function(e) { message("  H3K27ac VPA not loaded: ", e$message); GRanges() })

message("  H3K27ac CTRL peaks: ", length(h3k27_ctrl))
message("  H3K27ac VPA peaks:  ", length(h3k27_vpa))

# DMRs — all four pairwise contrasts
CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL", label="ASO alone",      col="#1F3A5F"),
  list(name="Scramble_VPA_vs_Scramble_CTRL", label="VPA alone",  col="#F0A500"),
  list(name="ASO_VPA_vs_Scramble_VPA",   label="ASO in VPA",     col="#C0392B"),
  list(name="ASO_VPA_vs_ASO_CTRL",       label="VPA in ASO",     col="#8E44AD")
)
dmr_list <- lapply(CONTRASTS, function(ct) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) return(GRanges())
  subsetByOverlaps(readRDS(rds), REGION)
})
names(dmr_list) <- sapply(CONTRASTS, function(x) x$label)
for (nm in names(dmr_list))
  message("  DMRs ", nm, ": ", length(dmr_list[[nm]]))


message("Loading methylation cache...")
meth_cache <- readRDS("results/dmr/meth_pooled_cache.rds")
meth_region <- lapply(meth_cache, function(m) subsetByOverlaps(m, REGION))


make_igv_plot <- function(outfile, title) {
  pdf(outfile, width=16, height=20, bg="white")

  # Track layout — r0/r1 values from bottom (0) to top (1)
  # Track 1 (top): gene model
  # Track 2: CpG islands
  # Track 3: H9 enhancers
  # Track 4: ENCODE cCREs
  # Track 5: H3K27ac CTRL
  # Track 6: H3K27ac VPA
  # Track 7-10: DMR tracks (one per contrast)
  # Track 11: methylation lines (4 conditions)

  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin     <- 0.12
  pp$topmargin      <- 50
  pp$bottommargin   <- 30
  pp$ideogramheight <- 8
  pp$data1height    <- 600

  kp <- plotKaryotype(
    genome      = "hg38",
    zoom        = REGION,
    plot.type   = 1,
    main        = title,
    cex         = 1.0,
    plot.params = pp
  )

  kpAddBaseNumbers(kp, tick.dist=50000, minor.tick.dist=10000, cex=0.7)

  
  # CpG islands
  if (length(cgi) > 0)
    kpPlotRegions(kp, data=cgi, col=CGI_COL, border=NA, r0=0.93, r1=0.98)
  kpAddLabels(kp, "CpG islands", r0=0.93, r1=0.98, cex=0.7, col=CGI_COL)

  # H9 enhancers
  if (length(h9_enh) > 0)
    kpPlotRegions(kp, data=h9_enh, col=ENH_H9_COL, border=NA, r0=0.86, r1=0.91)
  kpAddLabels(kp, "H9 ESC enh", r0=0.86, r1=0.91, cex=0.7, col=ENH_H9_COL)

  # ENCODE cCREs
  if (length(ccre) > 0)
    kpPlotRegions(kp, data=ccre, col=ENH_CCRE_COL, border=NA, r0=0.79, r1=0.84)
  kpAddLabels(kp, "ENCODE cCRE", r0=0.79, r1=0.84, cex=0.7, col=ENH_CCRE_COL)

  # H3K27ac CTRL
  if (length(h3k27_ctrl) > 0)
    kpPlotRegions(kp, data=h3k27_ctrl, col=H3K27_CTRL, border=NA, r0=0.72, r1=0.77)
  kpAddLabels(kp, "H3K27ac CTRL", r0=0.72, r1=0.77, cex=0.7, col=H3K27_CTRL)

  # H3K27ac VPA
  if (length(h3k27_vpa) > 0)
    kpPlotRegions(kp, data=h3k27_vpa, col=H3K27_VPA, border=NA, r0=0.65, r1=0.70)
  kpAddLabels(kp, "H3K27ac VPA", r0=0.65, r1=0.70, cex=0.7, col=H3K27_VPA)

  
  dmr_r0 <- c(0.57, 0.49, 0.41, 0.33)
  dmr_r1 <- c(0.63, 0.55, 0.47, 0.39)
  for (i in seq_along(CONTRASTS)) {
    ct   <- CONTRASTS[[i]]
    dmrs <- dmr_list[[ct$label]]
    # expand DMR regions by 5kb for visibility at this scale
    if (length(dmrs) > 0) {
      dmrs_wide <- dmrs
      start(dmrs_wide) <- pmax(1, start(dmrs_wide) - 5000)
      end(dmrs_wide)   <- end(dmrs_wide) + 5000
      kpPlotRegions(kp, data=dmrs_wide, col=adjustcolor(ct$col, 0.7),
                    border=ct$col, r0=dmr_r0[i], r1=dmr_r1[i])
    }
    kpAddLabels(kp, ct$label, r0=dmr_r0[i], r1=dmr_r1[i],
                cex=0.7, col=ct$col)
  }

  
  library(DMRcaller)
  for (cond in names(COND_COLS)) {
    tryCatch({
      win <- if (width(REGION) > 100000) 5000 else 1000
      prof <- computeMethylationProfile(meth_cache[[cond]], REGION,
                                        windowSize=win, context="CG")
      df   <- as.data.frame(prof)
      df   <- df[!is.na(df$sumReadsM) & df$sumReadsN >= 3, ]
      if (nrow(df) == 0) return(NULL)
      df$meth <- df$sumReadsM / df$sumReadsN
      df$pos  <- (df$start + df$end) / 2
      kpLines(kp, chr=as.character(seqnames(REGION)[1]),
              x=df$pos, y=df$meth,
              col=COND_COLS[cond], lwd=2.0, r0=0.0, r1=0.30)
    }, error=function(e) message("  meth profile failed for ", cond, ": ", e$message))
  }
  kpAddLabels(kp, "CpG methylation", r0=0.0, r1=0.30, cex=0.7, col="grey30")
  kpAxis(kp, r0=0.0, r1=0.30, ymin=0, ymax=1, numticks=3, cex=0.6)

  # Legend
  legend("topright", legend=names(COND_COLS), col=COND_COLS,
         lwd=2, bty="n", cex=0.8)

  # SMN2 gene boundary lines
  kpAbline(kp, v=SMN2_START, col="grey40", lwd=1.5, lty=2)
  kpAbline(kp, v=SMN2_END,   col="grey40", lwd=1.5, lty=2)

  dev.off()
  message("Saved: ", basename(outfile))
}

make_igv_plot(
  file.path(OUT, "SMN2_extended_IGV_200kb.pdf"),
  "SMN2 extended locus (±200kb) — methylation + DMRs + enhancers"
)

# Also make a tighter version ±50kb
REGION_TIGHT <- GRanges("chr5", IRanges(SMN2_START - 50000, SMN2_END + 50000))
meth_region  <- lapply(meth_cache, function(m) subsetByOverlaps(m, REGION_TIGHT))
cgi_raw2 <- read.table("data/reference/cpg_islands_hg38.bed", header=FALSE, sep="\t", stringsAsFactors=FALSE)
cgi <- subsetByOverlaps(GRanges(seqnames=cgi_raw2[,1], ranges=IRanges(cgi_raw2[,2], cgi_raw2[,3])), REGION_TIGHT)
h9_raw2 <- read.table(gzfile("data/reference/H9_predicted_non_promoter_non_fragments.bed.gz"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
h9_enh <- subsetByOverlaps(GRanges(seqnames=h9_raw2$seqnames, ranges=IRanges(h9_raw2$start, h9_raw2$end)), REGION_TIGHT)
ccre_raw3 <- read.table("data/reference/encode_cCREs_hg38.bed", header=FALSE, sep="\t", stringsAsFactors=FALSE)
ccre <- subsetByOverlaps(GRanges(seqnames=ccre_raw3[,1], ranges=IRanges(ccre_raw3[,2], ccre_raw3[,3])), REGION_TIGHT)
REGION       <- REGION_TIGHT

make_igv_plot(
  file.path(OUT, "SMN2_extended_IGV_50kb.pdf"),
  "SMN2 extended locus (±50kb) — methylation + DMRs + enhancers"
)

message("All done.")
