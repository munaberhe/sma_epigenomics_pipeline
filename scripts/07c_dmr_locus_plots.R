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
  list(name="RNA45SN2",  chr="chr21", start=8204909,   end=8213208,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="RNA45SN2 promoter (top ASO hit, p=1.56e-48)",
       annotation="Promoter (<=1kb) | hyper | -41.2%",
       gene_start=8208909, gene_end=8209208, strand="+",
       exons=data.frame(label="DMR", start=8208909, end=8209208, is_target=FALSE)),
  list(name="MTA1-DT",   chr="chr14", start=105412884, end=105423739,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="MTA1-DT intron (top hypo hit, p=1.72e-19)",
       annotation="Intron | hypo | +30.8%",
       gene_start=105414884, gene_end=105421739, strand="-",
       exons=data.frame(
         label=c("E1","E2","E3","E4","E5","E6"),
         start=c(105416884,105417581,105417833,105418682,105418939,105419022),
         end=c(105418309,105418309,105418312,105418816,105419080,105419739),
         is_target=rep(FALSE,6))),
  list(name="TRPV2",     chr="chr17", start=16403754,  end=16410053,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="TRPV2 (p=4.44e-22, axonogenesis GO)",
       annotation="Distal Intergenic | hyper | -33.5%",
       gene_start=16183742, gene_end=16268367, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="GLRA4",     chr="chrX",  start=103829053, end=103835352,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="GLRA4 chrX (p=5.86e-15, glycine receptor, neural)",
       annotation="Promoter (<=1kb) | hypo | chrX hotspot",
       gene_start=103829053, gene_end=103835352, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="NRXN1",     chr="chr2",  start=49914401,  end=50627200,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="NRXN1 (8 neural GO terms, synaptic)",
       annotation="Intron | 8 neural GO terms",
       gene_start=49914401, gene_end=50627200, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="PHACTR3",   chr="chr20", start=59601909,  end=59608208,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="PHACTR3 (p=1.02e-29, top ASO-on-VPA hypo)",
       annotation="Promoter | hypo | ASO_VPA vs Scramble_VPA",
       gene_start=59604540, gene_end=59847711, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="SOX5",      chr="chr12", start=24037319,  end=24043918,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="SOX5 (p=1.02e-13, neural TF, motor neuron)",
       annotation="Intron | hypo | motor neuron TF",
       gene_start=23533941, gene_end=24090418, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="SEMA3C",    chr="chr7",  start=80810801,  end=80817400,
       contrast="ASO_VPA_vs_Scramble_CTRL",
       label="SEMA3C (p=3.60e-69, axon guidance, neural)",
       annotation="Intron | hypo | axon guidance",
       gene_start=80765501, gene_end=80815949, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="HOXC8",     chr="chr12", start=54001919,  end=54008218,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="HOXC8 (p=3.41e-15, homeobox, motor neuron)",
       annotation="Promoter | hyper | motor neuron homeobox",
       gene_start=54004919, gene_end=54008985, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="SMN2",      chr="chr5",  start=70044638,  end=70083522,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="SMN2 therapeutic target (0 DMRs — null result)",
       annotation="Null result — no DMRs at therapeutic target",
       gene_start=70049638, gene_end=70078522, strand="-",
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
  # annotation label
  if (!is.null(locus$annotation)) {
    mtext(locus$annotation, side=3, line=0.2, cex=0.7,
          col="grey40", adj=1)
  }
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

# individual PDFs — one per locus using its specific contrast
message("\nindividual locus plots...")
for (locus in LOCI) {
  ct_name <- locus$contrast
  ct <- CONTRASTS[[which(sapply(CONTRASTS, function(x) x$name == ct_name))]]
  if (is.null(ct)) { message("  contrast not found: ", ct_name); next }
  {
    locus_gr   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
    dmrs_locus <- subsetByOverlaps(dmr_results[[ct$name]], locus_gr)
    message("  ", locus$name, " (", ct$name, "): ", length(dmrs_locus), " DMRs in region")
    out_path   <- file.path(OUT_DIR, paste0(locus$name,"_annotated.pdf"))
    cairo_pdf(out_path, width=11, height=6, bg="white")
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

# combined PDFs — locus-specific contrast only
message("\ncombined per-locus PDFs...")
for (locus in LOCI) {
  ct_name <- locus$contrast
  ct <- CONTRASTS[[which(sapply(CONTRASTS, function(x) x$name == ct_name))]]
  if (is.null(ct)) next
  locus_gr   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  dmrs_locus <- subsetByOverlaps(dmr_results[[ct$name]], locus_gr)
  out_path   <- file.path(OUT_DIR, paste0(locus$name, "_", ct$name, ".pdf"))
  cairo_pdf(out_path, width=11, height=6, bg="white")
  par(bg="white", col.axis="black", col.lab="black",
      col.main="black", fg="black", mar=c(5,4,4,2)+0.1, cex=0.85)
  tryCatch(
    plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
             ct, locus, dmrs_locus),
    error=function(e) { plot.new(); text(0.5,0.5,"failed",col="red") }
  )
  tryCatch(dev.off(), error=function(e) {
    message("  WARNING: dev.off failed for ", basename(out_path), ": ", e$message)
  })
  message("  saved: ", basename(out_path))
}
message("\ndone. outputs in: ", OUT_DIR)
