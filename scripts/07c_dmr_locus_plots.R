#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

# Annotated locus overlay plots for key DMR hits.
# Uses pooled methylation cache from 07b_dmr_plots.R.
# Exon coordinates from Ensembl GRCh38.109.
# SMN2 is handled separately by 05_smn2_locus_final.R.

DMR_DIR    <- "results/dmr"
OUT_DIR    <- "results/dmr/plots/annotated"
METH_CACHE <- file.path(DMR_DIR, "meth_pooled_cache.rds")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

WIN_SIZE <- 300

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond_a="ASO_CTRL", cond_b="Scramble_CTRL",
       label="ASO effect (CTRL background)"),
  list(name="ASO_VPA_vs_Scramble_CTRL",
       cond_a="ASO_VPA", cond_b="Scramble_CTRL",
       label="Combination vs baseline"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       cond_a="ASO_VPA", cond_b="ASO_CTRL",
       label="VPA effect on ASO background"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       cond_a="ASO_VPA", cond_b="Scramble_VPA",
       label="ASO effect on VPA background"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       cond_a="Scramble_VPA", cond_b="Scramble_CTRL",
       label="VPA effect (HDAC inhibitor)")
)

LOCI <- list(
  list(name="RNA45SN2", chr="chr21",
       start=8206909, end=8211208, strand="+",
       label="RNA45SN2 promoter (top ASO hit, p=1.56e-48)",
       gene_start=8208909, gene_end=8209208,
       exons=data.frame(label="DMR", start=8208909, end=8209208, is_target=TRUE)),
  list(name="MTA1-DT", chr="chr14",
       start=105414884, end=105421739, strand="-",
       label="MTA1-DT (top ASO hypo hit, intron, p=1.72e-19)",
       gene_start=105416884, gene_end=105419739,
       exons=data.frame(
         label=c("E1","E2","E3","E4","E5","E6"),
         start=c(105416884,105417581,105417833,105418682,105418939,105419022),
         end=c(105418309,105418309,105418312,105418816,105419080,105419739),
         is_target=rep(FALSE,6))),
  list(name="MYO1D", chr="chr17",
       start=32800000, end=32850000, strand="-",
       label="MYO1D (97 CpGs, p=7.83e-11, intronic)",
       gene_start=32492522, gene_end=32877177,
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0)))
)

# build GFF for plotLocalMethylationProfile gene track
build_gff <- function(locus) {
  rows <- list()
  rows[[1]] <- data.frame(chr=locus$chr, start=locus$gene_start,
    end=locus$gene_end, strand=locus$strand, type="gene", name=locus$name)
  for (i in seq_len(nrow(locus$exons)))
    rows[[length(rows)+1]] <- data.frame(chr=locus$chr,
      start=locus$exons$start[i], end=locus$exons$end[i],
      strand=locus$strand, type="exon",
      name=paste0(locus$name, "_", locus$exons$label[i]))
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          strand=df$strand, type=df$type, name=df$name)
}

plot_one <- function(meth_a, meth_b, ct, locus, dmrs=NULL) {
  region   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  gff      <- build_gff(locus)
  dmrs_arg <- if (!is.null(dmrs) && length(dmrs) > 0) list("DMRs"=dmrs) else NULL
  plotLocalMethylationProfile(
    methylationData1 = meth_a,
    methylationData2 = meth_b,
    region           = region,
    DMRs             = dmrs_arg,
    conditionsNames  = c(ct$cond_a, ct$cond_b),
    gff              = gff,
    windowSize       = WIN_SIZE,
    context          = "CG",
    main             = sprintf("%s: %s vs %s (%s)\n%s",
                               locus$name, ct$cond_a, ct$cond_b,
                               ct$label, locus$label),
    plotMeanLines = TRUE, plotPoints = TRUE
  )
  # exon labels — staggered to avoid overlap on dense genes
  ex <- locus$exons
  for (i in seq_len(nrow(ex))) {
    mtext(ex$label[i], side=1, at=(ex$start[i]+ex$end[i])/2,
          line=1.2, cex=0.65,
          col=if(ex$is_target[i]) "red" else "black",
          font=if(ex$is_target[i]) 2 else 1)
  }
  usr <- par("usr")
  text(usr[1], usr[3]+(usr[4]-usr[3])*0.05,
       labels=locus$name, cex=0.8, font=2, adj=c(0,0.5))
}

message("loading cached pooled data...")
if (!file.exists(METH_CACHE)) stop("cache not found: ", METH_CACHE)
meth_pooled <- readRDS(METH_CACHE)

message("loading DMR results...")
dmr_results <- list()
for (ct in CONTRASTS) {
  rds_path <- file.path(DMR_DIR, paste0("dmr_", ct$name, ".rds"))
  if (!file.exists(rds_path)) { message("  missing: ", ct$name); next }
  dmr_results[[ct$name]] <- readRDS(rds_path)
  message("  ", ct$name, ": ", length(dmr_results[[ct$name]]), " DMRs")
}

# individual PDFs — one per locus per contrast
message("\nindividual locus plots...")
for (ct in CONTRASTS) {
  for (locus in LOCI) {
    locus_gr   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
    dmrs_locus <- subsetByOverlaps(dmr_results[[ct$name]], locus_gr)
    out_path   <- file.path(OUT_DIR, paste0(ct$name,"_",locus$name,"_annotated.pdf"))
    pdf(out_path, width=11, height=6, bg="white")
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
    dev.off()
    message("  saved: ", basename(out_path))
  }
}

# combined PDFs — all contrasts for each locus
message("\ncombined per-locus PDFs...")
for (locus in LOCI) {
  locus_gr <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  out_path <- file.path(OUT_DIR, paste0(locus$name, "_all_contrasts.pdf"))
  pdf(out_path, width=16, height=14, bg="white")
  par(bg="white", col.axis="black", col.lab="black",
      col.main="black", fg="black", mfrow=c(3,1),
      mar=c(5,4,4,2)+0.1, cex=0.85)
  for (ct in CONTRASTS) {
    dmrs_locus <- subsetByOverlaps(dmr_results[[ct$name]], locus_gr)
    tryCatch(
      plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
               ct, locus, dmrs_locus),
      error=function(e) { plot.new(); text(0.5,0.5,"failed",col="red") }
    )
  }
  dev.off()
  message("  saved: ", basename(out_path))
}
message("\ndone. outputs in: ", OUT_DIR)
