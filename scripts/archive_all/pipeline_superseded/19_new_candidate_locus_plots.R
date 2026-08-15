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

# ASO context-dependent (key contrast: ASO_VPA_vs_Scramble_VPA)
# VPA context-dependent (key contrast: ASO_VPA_vs_ASO_CTRL)
LOCI <- list(
  list(name="IRF8",    chr="chr16", start=85865935, end=85966234,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="IRF8 — ASO context-dependent (transcription factor, chr16)",
       annotation="Promoter (1-2kb) | score=9 | ASO context-dependent",
       gene_start=85865935, gene_end=85966234, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="USP27X",  chr="chrX",  start=49832753, end=49933052,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="USP27X — ASO context-dependent (deubiquitinase, chrX)",
       annotation="Downstream | score=7 | ASO context-dependent",
       gene_start=49832753, gene_end=49933052, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="USP7",    chr="chr16", start=8840635,  end=8940934,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="USP7 — ASO context-dependent (deubiquitinase/chromatin, chr16)",
       annotation="Distal intergenic | score=6 | ASO context-dependent",
       gene_start=8840635,  gene_end=8940934,  strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="KDM1A",   chr="chr1",  start=23035569, end=23135868,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="KDM1A — ASO context-dependent (histone demethylase LSD1, chr1)",
       annotation="3' UTR | score=4 | histone demethylase",
       gene_start=23035569, gene_end=23135868, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="PAX5_ASO", chr="chr9", start=36992067, end=37092366,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="PAX5 — ASO context-dependent (transcription factor, chr9)",
       annotation="Intron | score=6 | ASO context-dependent",
       gene_start=36992067, gene_end=37092366, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="PAX5_VPA", chr="chr9", start=36992067, end=37092366,
       contrast="ASO_VPA_vs_ASO_CTRL",
       label="PAX5 — VPA context-dependent (transcription factor, chr9)",
       annotation="Intron | score=7 | VPA context-dependent",
       gene_start=36992067, gene_end=37092366, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="CAMK2A",  chr="chr5",  start=150164623, end=150264922,
       contrast="ASO_VPA_vs_ASO_CTRL",
       label="CAMK2A — VPA context-dependent (calcium kinase, motor neuron, chr5)",
       annotation="Intron | score=7 | synaptic/motor neuron",
       gene_start=150164623, gene_end=150264922, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="EPHB1",   chr="chr3",  start=135041420, end=135141719,
       contrast="ASO_VPA_vs_ASO_CTRL",
       label="EPHB1 — VPA context-dependent (axon guidance receptor, chr3)",
       annotation="Intron | score=7 | axon guidance",
       gene_start=135041420, gene_end=135141719, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="ZDHHC22", chr="chr14", start=77092523, end=77192822,
       contrast="ASO_VPA_vs_ASO_CTRL",
       label="ZDHHC22 — VPA context-dependent (palmitoyl transferase, chr14)",
       annotation="Promoter | score=6 | VPA context-dependent",
       gene_start=77092523, gene_end=77192822, strand="+",
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
  region <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  gff    <- build_gff(locus)
  meth_a <- meth_a[meth_a$readsN >= 10]
  meth_b <- meth_b[meth_b$readsN >= 10]

  plotLocalMethylationProfile(
    methylationData1 = meth_a,
    methylationData2 = meth_b,
    region           = region,
    DMRs             = NULL,
    conditionsNames  = c(ct$cond_a, ct$cond_b),
    gff              = gff,
    windowSize       = WIN_SIZE,
    context          = "CG",
    col              = c(COLS[ct$cond_a], COLS[ct$cond_b]),
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

message("\nPlotting new candidate loci...")
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
message("\nDone.")
