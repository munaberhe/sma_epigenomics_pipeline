#!/usr/bin/env Rscript
# 07j_karyoploter_bigwig.R
# ENCODE-style multi-track locus plots using kpPlotBigWig + autotrack,
# exactly following the karyoploteR ENCODE epigenetics tutorial pattern.
# Tracks (bottom to top):
#   Genes       -- kpPlotGenes via TxDb or manual segments
#   H9 enhancers -- kpPlotRegions (gold rectangles)
#   DMRs        -- kpPlotRegions (blue/red rectangles)
#   Methylation -- kpPlotBigWig per condition (cadetblue2 area tracks)
# Muna Berhe -- bt25018 -- QMUL MSc Bioinformatics

.libPaths('/data/home/bt25018/R/library')
suppressPackageStartupMessages({
  library(karyoploteR)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
})
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

BW_DIR  <- 'results/bigwigs'
DMR_DIR <- 'results/dmr'
H9_BED  <- 'data/ref/H9_predicted_non_promoter_non_fragments_bed.gz'
OUT_DIR <- 'results/dmr/plots/karyoploter'
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# Colours matching your deck
COND_COLS <- c(
  ASO_VPA       = '#F59E0B',
  ASO_CTRL      = '#02C39A',
  Scramble_CTRL = '#065A82',
  Scramble_VPA  = '#1C7293'
)
DMR_HYPO_COL  <- '#3182BD'
DMR_HYPER_COL <- '#E6550D'
ENH_COL       <- '#D4A017'

# Plot params -- matching tutorial style
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin    <- 0.15
pp$topmargin     <- 15
pp$bottommargin  <- 15
pp$ideogramheight <- 5
pp$data1inmargin  <- 10
pp$data1outmargin <- 0

# Load H9 enhancers once
message('Loading H9 enhancers...')
h9_dt <- fread(H9_BED, header=TRUE)
H9_ENH <- GRanges(seqnames=h9_dt$seqnames,
                  ranges=IRanges(h9_dt$start, h9_dt$end))
rm(h9_dt); gc()
message('  ', length(H9_ENH), ' enhancer regions loaded')

# Helper: load DMRs
load_dmrs <- function(contrast_name, chr, start, end) {
  rds <- file.path(DMR_DIR, paste0(contrast_name, '_DMRs.rds'))
  if (!file.exists(rds))
    rds <- file.path(DMR_DIR, paste0('dmr_', contrast_name, '.rds'))
  if (!file.exists(rds)) return(NULL)
  gr_all <- readRDS(rds)
  region <- GRanges(seqnames=chr, ranges=IRanges(start, end))
  gr <- subsetByOverlaps(gr_all, region)
  if (length(gr) == 0) return(NULL)
  if ('methylationDiff' %in% names(mcols(gr))) {
    mcols(gr)$direction <- ifelse(mcols(gr)$methylationDiff > 0, 'hypo', 'hyper')
  } else {
    mcols(gr)$direction <- 'hypo'
  }
  gr
}

# Helper: manual gene segments (no TxDb needed)
plot_genes <- function(kp, genes, r0, r1) {
  for (gene in genes) {
    gr <- GRanges(seqnames=gene$chr,
                  ranges=IRanges(gene$start, gene$end),
                  strand=gene$strand)
    kpSegments(kp, data=gr, y0=0.5, y1=0.5,
               r0=r0, r1=r1, col='grey25', lwd=2)
    mid <- (gene$start + gene$end) / 2
    kpText(kp, chr=gene$chr, x=mid, y=0.15,
           r0=r0, r1=r1,
           labels=gene$name, cex=1.2, col='grey10', font=2)
  }
}

# ---------------------------------------------------------------------------
# LOCI
# ---------------------------------------------------------------------------
LOCI <- list(

  list(
    name         = 'CACNG_cluster',
    chr          = 'chr17',
    start        = 66200000,
    end          = 67200000,
    title        = 'CACNG cluster -- chr17:66.2-67.2 Mb -- combination-exclusive hotspot',
    dmr_contrast = 'ASO_VPA_vs_Scramble_CTRL',
    conditions   = c('ASO_VPA', 'Scramble_CTRL'),
    genes = list(
      list(name='PRKCA',  chr='chr17', start=66236238, end=66595592, strand='-'),
      list(name='CACNG5', chr='chr17', start=66893661, end=66960714, strand='+'),
      list(name='CACNG4', chr='chr17', start=66979870, end=67046541, strand='+'),
      list(name='CACNG1', chr='chr17', start=67061729, end=67124036, strand='+')
    )
  ),

  list(
    name         = 'SEMA3C',
    chr          = 'chr7',
    start        = 80750000,
    end          = 80870000,
    title        = 'SEMA3C -- chr7:80.75-80.87 Mb -- axon guidance -- p=3.60e-69',
    dmr_contrast = 'ASO_VPA_vs_Scramble_CTRL',
    conditions   = c('ASO_VPA', 'Scramble_CTRL'),
    genes = list(
      list(name='SEMA3C', chr='chr7', start=80626327, end=81135324, strand='+')
    )
  ),

  list(
    name         = 'DNMBP',
    chr          = 'chr10',
    start        = 99900000,
    end          = 99930000,
    title        = 'DNMBP -- chr10:99.9-99.93 Mb -- promoter hypo -- p=4.57e-66',
    dmr_contrast = 'ASO_CTRL_vs_Scramble_CTRL',
    conditions   = c('ASO_CTRL', 'Scramble_CTRL'),
    genes = list(
      list(name='DNMBP', chr='chr10', start=99738672, end=100127694, strand='+')
    )
  ),

  list(
    name         = 'LINC00391_LMO7',
    chr          = 'chr13',
    start        = 75620000,
    end          = 94730000,
    title        = 'LMO7 / LINC00391 -- chr13 hotspot',
    dmr_contrast = 'ASO_CTRL_vs_Scramble_CTRL',
    conditions   = c('ASO_CTRL', 'Scramble_CTRL'),
    genes = list(
      list(name='LMO7',      chr='chr13', start=75528637, end=75781527, strand='+'),
      list(name='LINC00391', chr='chr13', start=94690000, end=94730000, strand='+')
    )
  ),

  list(
    name         = 'CHRNB3',
    chr          = 'chr8',
    start        = 42700000,
    end          = 42712000,
    title        = 'CHRNB3 -- chr8:42.70-42.71 Mb -- exon hypo -- p=8.18e-04',
    dmr_contrast = 'ASO_CTRL_vs_Scramble_CTRL',
    conditions   = c('ASO_CTRL', 'Scramble_CTRL'),
    genes = list(
      list(name='CHRNB3', chr='chr8', start=42607236, end=42830847, strand='+')
    )
  ),

  list(
    name         = 'HOXC8',
    chr          = 'chr12',
    start        = 54000000,
    end          = 54010000,
    title        = 'HOXC8 -- chr12:54.00-54.01 Mb -- promoter hyper -- p=3.41e-15',
    dmr_contrast = 'ASO_VPA_vs_Scramble_VPA',
    conditions   = c('ASO_VPA', 'Scramble_VPA'),
    genes = list(
      list(name='HOXC8', chr='chr12', start=53974956, end=54002785, strand='-')
    )
  )
)

# ---------------------------------------------------------------------------
# Plot function -- ENCODE tutorial pattern
# ---------------------------------------------------------------------------
plot_locus <- function(locus) {
  region_gr <- toGRanges(paste0(locus$chr, ':', locus$start, '-', locus$end))
  n_cond    <- length(locus$conditions)

  # Check bigwigs exist
  bw_files <- sapply(locus$conditions, function(c)
    file.path(BW_DIR, paste0(c, '_methylation.bw')))
  bw_ok <- file.exists(bw_files)
  if (!any(bw_ok)) {
    message('  SKIP -- no bigwig files found. Run 07i_export_bigwigs.R first.')
    return(invisible(NULL))
  }

  outfile <- file.path(OUT_DIR, paste0(locus$name, '_karyoploter.pdf'))
  pdf(outfile, width=12, height=5 + n_cond * 1.5)

  kp <- plotKaryotype(
    zoom        = region_gr,
    genome      = 'hg38',
    plot.type   = 1,
    cex         = 2,
    plot.params = pp
  )
  kpAddMainTitle(kp, locus$title, cex=1.8)
  kpAddBaseNumbers(kp, tick.dist=50000, minor.tick.dist=10000,
                   add.units=TRUE, cex=1.2)

  # Track r0/r1 layout (bottom to top):
  #   0.00 - 0.10  genes
  #   0.12 - 0.17  H9 enhancers
  #   0.19 - 0.25  DMRs
  #   0.28 - 1.00  methylation bigwig tracks (autotrack)
  gene_r0 <- 0.00; gene_r1 <- 0.10
  enh_r0  <- 0.12; enh_r1  <- 0.17
  dmr_r0  <- 0.19; dmr_r1  <- 0.25
  meth_r0 <- 0.28; meth_r1 <- 1.00

  # ---- Gene track ----
  plot_genes(kp, locus$genes, gene_r0, gene_r1)
  kpAddLabels(kp, labels='Genes', r0=gene_r0, r1=gene_r1,
              cex=1.4, label.margin=0.035)

  # ---- H9 enhancers ----
  enh_in <- subsetByOverlaps(H9_ENH, region_gr)
  if (length(enh_in) > 0) {
    kpPlotRegions(kp, enh_in,
                  col=adjustcolor(ENH_COL, alpha.f=0.7),
                  r0=enh_r0, r1=enh_r1, border=NA)
    message('  H9 enhancers in window: ', length(enh_in))
  }
  kpAddLabels(kp, labels='H9 enh', r0=enh_r0, r1=enh_r1,
              cex=1.4, col=ENH_COL, label.margin=0.035)

  # ---- DMR rectangles ----
  dmr_gr <- load_dmrs(locus$dmr_contrast, locus$chr, locus$start, locus$end)
  if (!is.null(dmr_gr) && length(dmr_gr) > 0) {
    dmr_cols <- ifelse(mcols(dmr_gr)$direction == 'hypo',
                       DMR_HYPO_COL, DMR_HYPER_COL)
    kpPlotRegions(kp, dmr_gr, col=dmr_cols,
                  r0=dmr_r0, r1=dmr_r1, border=NA)
    message('  DMRs plotted: ', length(dmr_gr))
  }
  kpAddLabels(kp, labels='DMRs', r0=dmr_r0, r1=dmr_r1,
              cex=1.4, col='grey30', label.margin=0.035)

  # ---- Methylation bigwig tracks (autotrack -- tutorial pattern) ----
  valid_conds <- locus$conditions[bw_ok]
  for (i in seq_along(valid_conds)) {
    cond <- valid_conds[i]
    col  <- COND_COLS[cond]
    bw   <- bw_files[cond]
    at   <- autotrack(i, length(valid_conds),
                      r0=meth_r0, r1=meth_r1, margin=0.1)
    kp <- kpPlotBigWig(kp, data=bw,
                       ymax='visible.region',
                       r0=at$r0, r1=at$r1,
                       col=adjustcolor(col, alpha.f=0.7),
                       border=col)
    computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax * 10) / 10
    kpAxis(kp, ymin=0, ymax=computed.ymax,
           tick.pos=c(0, computed.ymax),
           numticks=2, r0=at$r0, r1=at$r1, cex=1.2)
    kpAddLabels(kp, labels=cond, r0=at$r0, r1=at$r1,
                cex=1.4, col=col, label.margin=0.035)
  }

  # Legend
  legend('topright',
         legend=c('Hypo DMR', 'Hyper DMR', 'H9 enhancer'),
         fill=c(DMR_HYPO_COL, DMR_HYPER_COL,
                adjustcolor(ENH_COL, alpha.f=0.7)),
         bty='n', cex=1.2)

  dev.off()
  message('  Saved: ', basename(outfile))
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
message('\nPlotting loci...')
for (locus in LOCI) {
  message('\nLocus: ', locus$name)
  tryCatch(
    plot_locus(locus),
    error = function(e) message('  ERROR: ', e$message)
  )
}
message('\nDone. Outputs in: ', OUT_DIR)
