#!/usr/bin/env Rscript
# 07e_karyoploter_tracks.R
# Multi-track karyoploteR locus plots replicating the ENCODE epigenetics style:
#   Track 1 — CpG methylation proportion per condition (area plots, one per condition)
#   Track 2 — DMR rectangles coloured by direction (hypo = blue, hyper = red)
#   Track 3 — Gene / exon model
# Uses the same meth_pooled_cache.rds and DMR .rds files as 07c/07d.
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(karyoploteR)
})
.libPaths(c("~/R/library", .libPaths()))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DMR_DIR    <- "results/dmr"
OUT_DIR    <- "results/dmr/plots/karyoploter"
METH_CACHE <- file.path(DMR_DIR, "meth_pooled_cache.rds")

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ---------------------------------------------------------------------------
# Colours — match your existing deck palette
# ---------------------------------------------------------------------------
COND_COLS <- c(
  ASO_CTRL      = "#02C39A",   # mint
  Scramble_CTRL = "#065A82",   # deep blue
  ASO_VPA       = "#F59E0B",   # amber
  Scramble_VPA  = "#1C7293"    # teal
)
DMR_HYPO_COL  <- "#3182BD"   # blue
DMR_HYPER_COL <- "#E6550D"   # red-orange

# ---------------------------------------------------------------------------
# Loci to plot
# Each locus specifies which contrast drives the DMR annotation and which
# conditions to show in the methylation tracks.
# ---------------------------------------------------------------------------
LOCI <- list(

  list(
    name       = "CACNG_cluster",
    chr        = "chr17",
    start      = 66200000,
    end        = 67200000,
    dmr_contrast = "ASO_VPA_vs_Scramble_CTRL",
    conditions   = c("ASO_VPA", "Scramble_CTRL"),
    label      = "CACNG cluster · chr17:66.2–67.2 Mb · combination-exclusive hotspot",
    genes = list(
      list(name="PRKCA",  start=66236238, end=66595592, strand="-"),
      list(name="CACNG5", start=66893661, end=66960714, strand="+"),
      list(name="CACNG4", start=66979870, end=67046541, strand="+"),
      list(name="CACNG1", start=67061729, end=67124036, strand="+")
    ),
    dmr_highlights = list(
      list(start=66311454, end=66311753, gene="PRKCA",  direction="hyper"),
      list(start=66911154, end=66911753, gene="CACNG5", direction="hypo"),
      list(start=66997554, end=66998153, gene="CACNG4", direction="hyper"),
      list(start=67079454, end=67079753, gene="CACNG1", direction="hyper")
    )
  ),

  list(
    name       = "SEMA3C",
    chr        = "chr7",
    start      = 80750000,
    end        = 80870000,
    dmr_contrast = "ASO_VPA_vs_Scramble_CTRL",
    conditions   = c("ASO_VPA", "Scramble_CTRL"),
    label      = "SEMA3C · chr7:80.75–80.87 Mb · axon guidance · p=3.60e-69",
    genes = list(
      list(name="SEMA3C", start=80626327, end=81135324, strand="+")
    ),
    dmr_highlights = list(
      list(start=80813801, end=80814400, gene="SEMA3C", direction="hypo")
    )
  ),

  list(
    name       = "DNMBP",
    chr        = "chr10",
    start      = 99900000,
    end        = 99930000,
    dmr_contrast = "ASO_CTRL_vs_Scramble_CTRL",
    conditions   = c("ASO_CTRL", "Scramble_CTRL"),
    label      = "DNMBP · chr10:99.9–99.93 Mb · promoter hypo · p=4.57e-66",
    genes = list(
      list(name="DNMBP", start=99738672, end=100127694, strand="+")
    ),
    dmr_highlights = list(
      list(start=99914000, end=99920000, gene="DNMBP", direction="hypo")
    )
  ),

  list(
    name       = "LINC00391_LMO7",
    chr        = "chr13",
    start      = 75620000,
    end        = 94730000,
    dmr_contrast = "ASO_CTRL_vs_Scramble_CTRL",
    conditions   = c("ASO_CTRL", "Scramble_CTRL"),
    label      = "LMO7 / LINC00391 · chr13 hotspot",
    genes = list(
      list(name="LMO7",      start=75528637, end=75781527, strand="+"),
      list(name="LINC00391", start=94690000, end=94730000, strand="+")
    ),
    dmr_highlights = list(
      list(start=75632000, end=75638000, gene="LMO7",      direction="hypo"),
      list(start=94707000, end=94713000, gene="LINC00391", direction="hypo")
    )
  ),

  list(
    name         = "CHRNB3",
    chr          = "chr8",
    start        = 42700000,
    end          = 42712000,
    dmr_contrast = "ASO_CTRL_vs_Scramble_CTRL",
    conditions   = c("ASO_CTRL", "Scramble_CTRL"),
    label        = "CHRNB3 chr8:42.70-42.71 Mb exon hypo p=8.18e-04",
    genes = list(
      list(name="CHRNB3", start=42607236, end=42830847, strand="+")
    ),
    dmr_highlights = list(
      list(start=42703000, end=42709000, gene="CHRNB3", direction="hypo")
    )
  ),

  list(
    name         = "HOXC8",
    chr          = "chr12",
    start        = 54000000,
    end          = 54010000,
    dmr_contrast = "ASO_VPA_vs_Scramble_VPA",
    conditions   = c("ASO_VPA", "Scramble_VPA"),
    label        = "HOXC8 chr12:54.00-54.01 Mb promoter hyper p=3.41e-15",
    genes = list(
      list(name="HOXC8", start=53974956, end=54002785, strand="-")
    ),
    dmr_highlights = list(
      list(start=54002000, end=54008000, gene="HOXC8", direction="hyper")
    )
  )

)

# ---------------------------------------------------------------------------
# Helper: extract per-CpG methylation from a pooled DMRcaller object as GRanges
# with a score column (methylation proportion)
# ---------------------------------------------------------------------------
meth_to_gr <- function(meth_obj, chr, start, end) {
  # DMRcaller pooled objects are GRanges with readsM and readsN columns
  region <- GRanges(seqnames=chr, ranges=IRanges(start, end))
  hits   <- subsetByOverlaps(meth_obj, region)
  if (length(hits) == 0) return(NULL)
  # Compute methylation proportion; guard against 0 coverage
  cov  <- mcols(hits)$readsN
  m    <- mcols(hits)$readsM
  prop <- ifelse(cov > 0, m / cov, NA_real_)
  keep <- !is.na(prop)
  hits <- hits[keep]
  mcols(hits)$score <- prop[keep]
  hits
}

# ---------------------------------------------------------------------------
# Helper: build DMR GRanges coloured by direction from the annotated CSV
# Falls back to the raw DMR .rds if the CSV is absent.
# ---------------------------------------------------------------------------
load_dmrs_for_locus <- function(contrast_name, chr, start, end) {
  csv_path <- file.path(DMR_DIR, "annotation",
                        paste0(contrast_name, "_annotated.csv"))
  rds_path <- file.path(DMR_DIR, paste0(contrast_name, "_DMRs.rds"))
  if (!file.exists(rds_path))
    rds_path <- file.path(DMR_DIR, paste0("dmr_", contrast_name, ".rds"))

  if (file.exists(csv_path)) {
    df  <- read.csv(csv_path, stringsAsFactors=FALSE)
    df  <- df[df$seqnames == chr &
              df$start    >= start &
              df$end      <= end, ]
    if (nrow(df) == 0) return(NULL)
    gr  <- GRanges(seqnames=df$seqnames,
                   ranges=IRanges(df$start, df$end))
    # regionType column: "gain" = hypomethylation in DMRcaller convention
    mcols(gr)$direction <- ifelse(df$regionType == "gain", "hypo", "hyper")
    mcols(gr)$gene      <- if ("SYMBOL" %in% names(df)) df$SYMBOL else ""
    return(gr)
  } else if (file.exists(rds_path)) {
    gr_all <- readRDS(rds_path)
    region <- GRanges(seqnames=chr, ranges=IRanges(start, end))
    gr     <- subsetByOverlaps(gr_all, region)
    if (length(gr) == 0) return(NULL)
    # Infer direction from methylationDiff column if present
    if ("methylationDiff" %in% names(mcols(gr))) {
      mcols(gr)$direction <- ifelse(mcols(gr)$methylationDiff > 0, "gain", "loss")
    } else {
      mcols(gr)$direction <- "hypo"
    }
    return(gr)
  }
  return(NULL)
}

# ---------------------------------------------------------------------------
# Plot one locus
# ---------------------------------------------------------------------------
plot_locus <- function(locus, meth_pooled) {

  region_gr <- toGRanges(paste0(locus$chr, ":",
                                locus$start, "-", locus$end))

  # How many methylation tracks?
  n_meth  <- length(locus$conditions)
  # Track layout (r0 to r1, bottom to top):
  #   0.00 – 0.10  gene model
  #   0.10 – 0.18  DMR bar
  #   0.18 – 0.20  spacer
  #   0.20 – 1.00  methylation tracks (equal split)
  gene_r0 <- 0.00;  gene_r1 <- 0.10
  dmr_r0  <- 0.10;  dmr_r1  <- 0.18
  meth_height <- (1.00 - 0.20) / n_meth
  meth_r0 <- function(i) 0.20 + (i - 1) * meth_height
  meth_r1 <- function(i) 0.20 +  i      * meth_height

  outfile <- file.path(OUT_DIR, paste0(locus$name, "_karyoploter.pdf"))
  pdf(outfile, width=11, height=5 + n_meth * 0.8)

  kp <- plotKaryotype(
    zoom         = region_gr,
    genome       = "hg38",
    plot.type    = 1,
    cex          = 1.2,
    main         = locus$label
  )

  # ---- Methylation area tracks (one per condition) ------------------------
  for (i in seq_along(locus$conditions)) {
    cond  <- locus$conditions[i]
    col   <- COND_COLS[cond]
    m_obj <- meth_pooled[[cond]]
    if (is.null(m_obj)) {
      message("  Condition not in cache: ", cond)
      next
    }
    cpg_gr <- meth_to_gr(m_obj, locus$chr, locus$start, locus$end)
    if (is.null(cpg_gr) || length(cpg_gr) == 0) {
      message("  No CpGs for ", cond, " in window")
      next
    }
    r0 <- meth_r0(i); r1 <- meth_r1(i)
    kpArea(kp, data=cpg_gr, y=mcols(cpg_gr)$score,
           r0=r0, r1=r1,
           col=adjustcolor(col, alpha.f=0.6),
           border=col,
           ymin=0, ymax=1)
    kpAxis(kp, r0=r0, r1=r1, ymin=0, ymax=1,
           tick.pos=c(0, 0.5, 1), cex=0.5, col="grey50")
    kpAddLabels(kp, labels=cond, r0=r0, r1=r1,
                cex=0.65, col=col, label.margin=0.03)
  }

  # ---- DMR rectangles ----------------------------------------------------
  dmr_gr <- load_dmrs_for_locus(locus$dmr_contrast,
                                 locus$chr, locus$start, locus$end)
  if (!is.null(dmr_gr) && length(dmr_gr) > 0) {
    dmr_cols <- ifelse(mcols(dmr_gr)$direction == "hypo",
                       DMR_HYPO_COL, DMR_HYPER_COL)
    kpRect(kp, data=dmr_gr,
           y0=0, y1=1,
           r0=dmr_r0, r1=dmr_r1,
           col=dmr_cols, border=NA)
  }
  # DMR track label
  kpAddLabels(kp, labels="DMRs", r0=dmr_r0, r1=dmr_r1,
              cex=0.6, col="grey30", label.margin=0.03)
  # DMR legend
  legend("topright",
         legend=c("Hypo DMR", "Hyper DMR"),
         fill=c(DMR_HYPO_COL, DMR_HYPER_COL),
         bty="n", cex=0.65)

  # ---- Gene track --------------------------------------------------------
  for (gene in locus$genes) {
    gene_gr <- GRanges(seqnames=locus$chr,
                       ranges=IRanges(gene$start, gene$end),
                       strand=gene$strand)
    # Gene body line
    kpSegments(kp, data=gene_gr,
               y0=0.5, y1=0.5,
               r0=gene_r0, r1=gene_r1,
               col="grey30", lwd=1.5)
    # Gene name label
    mid <- (gene$start + gene$end) / 2
    kpText(kp,
           chr=locus$chr, x=mid, y=0.15,
           r0=gene_r0, r1=gene_r1,
           labels=gene$name,
           cex=0.6, col="grey10")
  }
  kpAddLabels(kp, labels="Genes", r0=gene_r0, r1=gene_r1,
              cex=0.6, col="grey30", label.margin=0.03)

  dev.off()
  message("  Saved: ", basename(outfile))
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
message("Loading methylation cache...")
meth_pooled <- readRDS(METH_CACHE)

message("Plotting loci...")
for (locus in LOCI) {
  message("\nLocus: ", locus$name)
  tryCatch(
    plot_locus(locus, meth_pooled),
    error = function(e) message("  ERROR: ", e$message)
  )
}

message("\nDone. Output in: ", OUT_DIR)
