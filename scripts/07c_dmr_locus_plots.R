#!/usr/bin/env Rscript
# 07c_dmr_locus_plots.R
# Annotated locus overlay plots for key DMR hits
# Follows the same approach as 05_smn2_locus_final.R
# DMRs passed as named list per DMRcaller vignette section 3.15
# Exon coordinates hardcoded from Ensembl GRCh38.109

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
.libPaths(c("~/R/library", .libPaths()))

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

# Locus definitions with hardcoded gene body and exon coordinates.
# Exon boundaries from Ensembl GRCh38.109, deduplicated to unique start/end pairs.
# SMN2 excluded -- 0 DMRs at locus across all contrasts, dedicated masked script
# 05_smn2_locus_final.R produces better plots for that locus.
# RNA45SN2 has no exon annotation in Ensembl (rDNA locus) -- DMR window marked only.
# MTA1-DT exons are densely packed (min gap 127bp) so labels are staggered.
# MYO1D DMR falls in an intron -- no exons in the 32,800,000-32,850,000 window.
LOCI <- list(

  list(
    name       = "RNA45SN2",
    chr        = "chr21",
    start      = 8206909,
    end        = 8211208,
    strand     = "+",
    label      = "RNA45SN2 promoter (top ASO hit, p=1.56e-48, 109 CpGs, 41% drop)",
    gene_start = 8208909,
    gene_end   = 8209208,
    exons = data.frame(
      label     = "DMR",
      start     = 8208909,
      end       = 8209208,
      is_target = TRUE
    )
  ),

  list(
    name       = "MTA1-DT",
    chr        = "chr14",
    start      = 105414884,
    end        = 105421739,
    strand     = "-",
    label      = "MTA1-DT locus (top ASO hypo hit, intron, p=1.72e-19)",
    gene_start = 105416884,
    gene_end   = 105419739,
    exons = data.frame(
      label     = c("E1","E2","E3","E4","E5","E6"),
      start     = c(105416884,105417581,105417833,105418682,105418939,105419022),
      end       = c(105418309,105418309,105418312,105418816,105419080,105419739),
      is_target = c(F,F,F,F,F,F)
    )
  ),

  list(
    name       = "MYO1D",
    chr        = "chr17",
    start      = 32800000,
    end        = 32850000,
    strand     = "-",
    label      = "MYO1D locus (97 CpGs, p=7.83e-11) -- DMR at intron",
    gene_start = 32492522,
    gene_end   = 32877177,
    exons = data.frame(
      label     = character(0),
      start     = numeric(0),
      end       = numeric(0),
      is_target = logical(0)
    )
  )
)

# Build GFF from hardcoded coordinates for plotLocalMethylationProfile.
# Requires GRanges with type and name metadata columns per vignette section 3.15.
build_gff <- function(locus) {
  rows <- list()
  rows[[1]] <- data.frame(
    chr=locus$chr, start=locus$gene_start, end=locus$gene_end,
    strand=locus$strand, type="gene", name=locus$name
  )
  for (i in seq_len(nrow(locus$exons))) {
    rows[[length(rows)+1]] <- data.frame(
      chr=locus$chr, start=locus$exons$start[i], end=locus$exons$end[i],
      strand=locus$strand, type="exon",
      name=paste0(locus$name, "_", locus$exons$label[i])
    )
  }
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          strand=df$strand, type=df$type, name=df$name)
}

# Plot one locus for one contrast.
# DMRs passed as named list per vignette: list("DMRs"=GRanges).
# Exon labels added with mtext() -- alternating line heights to avoid overlap
# on densely packed genes like MTA1-DT (min inter-exon gap = 127bp).
plot_one <- function(meth_a, meth_b, ct, locus, dmrs=NULL) {
  region <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  gff    <- build_gff(locus)

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
    plotMeanLines    = TRUE,
    plotPoints       = TRUE
  )

  # Add exon labels below gene track -- alternating line heights (0.5 and 1.8)
  # to prevent overlap on densely packed genes
  ex <- locus$exons
  for (i in seq_len(nrow(ex))) {
    mid     <- (ex$start[i] + ex$end[i]) / 2
    is_tgt  <- ex$is_target[i]
    stagger <- if (i %% 2 == 0) 0.5 else 1.8
    mtext(ex$label[i], side=1, at=mid, line=stagger,
          cex=0.45, col=if(is_tgt) "red" else "black",
          font=if(is_tgt) 2 else 1)
  }

  # Gene name at bottom left
  usr <- par("usr")
  text(usr[1], usr[3] + (usr[4]-usr[3])*0.05,
       labels=locus$name, cex=0.8, font=2, col="black", adj=c(0,0.5))
}

# Load pooled methylation cache from 07b_dmr_plots.R
message("Loading cached pooled methylation data...")
if (!file.exists(METH_CACHE)) stop("Cache not found: ", METH_CACHE)
meth_pooled <- readRDS(METH_CACHE)
message("  Loaded: ", paste(names(meth_pooled), collapse=", "))

# Load DMR results so called DMRs can be shown as boxes on the plots
message("Loading DMR results...")
dmr_results <- list()
for (ct in CONTRASTS) {
  rds_path <- file.path(DMR_DIR, paste0("dmr_", ct$name, ".rds"))
  if (!file.exists(rds_path)) { message("  Missing: ", ct$name); next }
  dmr_results[[ct$name]] <- readRDS(rds_path)
  message("  ", ct$name, ": ", length(dmr_results[[ct$name]]), " DMRs")
}

# Individual PDFs -- one per locus per contrast
message("\nGenerating individual locus plots...")

for (ct in CONTRASTS) {
  meth_a <- meth_pooled[[ct$cond_a]]
  meth_b <- meth_pooled[[ct$cond_b]]

  for (locus in LOCI) {
    message("  ", locus$name, " / ", ct$name)

    locus_gr   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
    dmrs_locus <- subsetByOverlaps(dmr_results[[ct$name]], locus_gr)
    message("    DMRs in window: ", length(dmrs_locus))

    out_path <- file.path(OUT_DIR, paste0(ct$name, "_", locus$name, "_annotated.pdf"))

    pdf(out_path, width=11, height=6, bg="white")
    par(bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")
    par(mar=c(5, 4, 4, 2) + 0.1)

    tryCatch(
      plot_one(meth_a, meth_b, ct, locus, dmrs_locus),
      error = function(e) {
        message("    Error: ", conditionMessage(e))
        plot.new()
        text(0.5, 0.5,
             paste0("Plot failed\n", conditionMessage(e),
                    "\nLocus: ", locus$name, "\nContrast: ", ct$name),
             col="red", cex=0.9)
      }
    )

    dev.off()
    message("    Saved: ", basename(out_path))
  }
}

# Combined PDF -- all three contrasts stacked for each locus.
# Most useful for thesis -- shows how each treatment affects the same locus.
message("\nGenerating combined per-locus PDFs...")

for (locus in LOCI) {
  message("  Combined: ", locus$name)
  locus_gr <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  out_path <- file.path(OUT_DIR, paste0(locus$name, "_all_contrasts.pdf"))

  pdf(out_path, width=16, height=14, bg="white")
  par(bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")
  par(mfrow=c(3,1), mar=c(5, 4, 4, 2) + 0.1, cex=0.85)

  for (ct in CONTRASTS) {
    dmrs_locus <- subsetByOverlaps(dmr_results[[ct$name]], locus_gr)
    tryCatch(
      plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
               ct, locus, dmrs_locus),
      error = function(e) {
        message("    Error: ", conditionMessage(e))
        plot.new()
        text(0.5, 0.5, paste0("Failed: ", conditionMessage(e)), col="red")
      }
    )
  }

  dev.off()
  message("    Saved: ", basename(out_path))
}

message("\nDone. Outputs in: ", OUT_DIR)
plots <- list.files(OUT_DIR, pattern="\\.pdf$")
for (f in sort(plots)) message("  ", f)
