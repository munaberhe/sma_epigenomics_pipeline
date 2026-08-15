#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
source("scripts/pipeline/00_sma_palette.R")

OUT_DIR    <- "results/dmr/plots/annotated"
METH_CACHE <- "results/dmr/meth_pooled_cache.rds"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
WIN_SIZE <- 300

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond_a="Scramble_CTRL", cond_b="ASO_CTRL",
       label="ASO effect (CTRL background)"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       cond_a="Scramble_CTRL", cond_b="Scramble_VPA",
       label="VPA effect (HDAC inhibitor)"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       cond_a="Scramble_VPA", cond_b="ASO_VPA",
       label="ASO effect on VPA background"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       cond_a="ASO_CTRL", cond_b="ASO_VPA",
       label="VPA effect on ASO background")
)

# 50kb window around each gene
LOCI <- list(
  list(name="SEMA3C", chr="chr7",
       start=80760801, end=80867400,
       contrast="Scramble_VPA_vs_Scramble_CTRL",
       label="SEMA3C — VPA-driven example (axon guidance, chr7)",
       annotation="Intron | hypo | VPA-driven | axon guidance",
       gene_start=80765501, gene_end=80815949, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="SIGMAR1", chr="chr9",
       start=34589667, end=34692366,
       contrast="ASO_VPA_vs_ASO_CTRL",
       label="SIGMAR1 — ASO-driven example (sigma receptor, chr9)",
       annotation="Promoter | hyper | ASO-driven | sigma receptor",
       gene_start=34634739, gene_end=34637844, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0)))
)

build_gff <- function(locus) {
  rows <- list()
  rows[[1]] <- data.frame(chr=locus$chr, start=locus$gene_start,
                          end=locus$gene_end, strand=locus$strand,
                          type="gene", name=locus$name)
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          strand=df$strand, type=df$type, name=df$name)
}

COLS <- c(ASO_CTRL="#1B4F8A", ASO_VPA="#B2182B",
          Scramble_VPA="#F0A500", Scramble_CTRL="#6B7280")

plot_one <- function(meth_a, meth_b, ct, locus, dmrs=NULL) {
  region   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  gff      <- build_gff(locus)
  meth_a   <- meth_a[meth_a$readsN >= 10]
  meth_b   <- meth_b[meth_b$readsN >= 10]
  plotLocalMethylationProfile(
    methylationData1 = meth_a,
    methylationData2 = meth_b,
    region           = region,
    DMRs             = NULL,
    conditionsNames  = c(ct$cond_a, ct$cond_b),
    gff              = gff,
    windowSize       = WIN_SIZE,
    context          = "CG",
    col              = c(COLS[ct$cond_a], COLS[ct$cond_b],
                         "#1F3A5F", "#6B7280", "#E31A1C"),
    main             = sprintf("%s: %s vs %s (%s)\n%s",
                               locus$name, ct$cond_a, ct$cond_b,
                               ct$label, locus$label),
    plotMeanLines = TRUE, plotPoints = FALSE
  )
  # overdraw DMR positions as narrow red band at top of methylation panel
  if (!is.null(dmrs) && length(dmrs) > 0 && length(dmrs) <= 30) {
    # switch back to methylation panel (panel 1) before drawing
    tryCatch(par(mfg=c(1,1)), error=function(e) NULL)
    usr <- par("usr")
    meth_top <- usr[4] - (usr[4]-usr[3]) * 0.05
    meth_bot <- usr[4] - (usr[4]-usr[3]) * 0.12
    dmr_df <- as.data.frame(dmrs)
    pad <- 1000
    for (i in seq_len(nrow(dmr_df))) {
      rect(xleft   = dmr_df$start[i] - pad,
           xright  = dmr_df$end[i]   + pad,
           ybottom = meth_bot,
           ytop    = meth_top,
           col     = "#E31A1C",
           border  = "#E31A1C",
           lwd     = 0.5)
    }
  }
  if (!is.null(locus$annotation))
    mtext(locus$annotation, side=3, line=0.2, cex=0.7, col="grey40", adj=1)
}

message("Loading pooled methylation cache...")
meth_pooled <- readRDS(METH_CACHE)

message("Loading DMR results...")
dmr_results <- list()
for (ct in CONTRASTS) {
  rds <- file.path("results/dmr", paste0("dmr_", ct$name, ".rds"))
  if (file.exists(rds)) {
    dmr_results[[ct$name]] <- readRDS(rds)
    message("  ", ct$name, ": ", length(dmr_results[[ct$name]]), " DMRs")
  }
}

message("\nPlotting SEMA3C and SIGMAR1...")
for (locus in LOCI) {
  out_path <- file.path(OUT_DIR, paste0(locus$name, "_4contrasts.pdf"))
  cairo_pdf(out_path, width=11, height=6, bg="white", onefile=TRUE)
  for (ct in CONTRASTS) {
    if (is.null(dmr_results[[ct$name]])) next
    locus_gr   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
    dmrs_locus <- subsetByOverlaps(dmr_results[[ct$name]], locus_gr)
    par(bg="white", col.axis="black", col.lab="black",
        col.main="black", fg="black", mar=c(5,4,4,2)+0.1)
    tryCatch(
      plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
               ct, locus, dmrs_locus),
      error=function(e) {
        plot.new()
        text(0.5, 0.5, paste0("plot failed: ", conditionMessage(e)),
             col="red", cex=0.9)
      }
    )
  }
  dev.off()
  message("Saved: ", basename(out_path))
}
message("Done.")
