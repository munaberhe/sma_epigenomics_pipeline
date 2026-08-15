#!/usr/bin/env Rscript
# 21_intron67_enhancer_chr13_annotation.R
# Item 9: Check H3K27ac peaks specifically in SMN2 introns 6-7
# Item 10: Annotate chromosome 13 DMR hotspot genes

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ChIPseeker)
})
.libPaths(c('~/R/library', .libPaths()))
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

# ── ITEM 9: SMN2 introns 6 and 7 H3K27ac check ───────────────────────────────
message('=== ITEM 9: SMN2 intron 6-7 H3K27ac analysis ===')

# SMN2 exon coordinates (Alberto's convention, E7 = penultimate exon)
# Intron 6 = between E6 and E7
# Intron 7 = between E7 and E8
EXONS <- data.frame(
  exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
  start = c(70049638,70053107,70056229,70063044,70069090,
            70069235,70070641,70076521,70077019),
  end   = c(70050437,70053264,70056357,70063153,70069186,
            70069330,70070751,70076574,70077592)
)

# Intron boundaries
intron6_start <- EXONS$end[EXONS$exon=="E6"]   # after E6
intron6_end   <- EXONS$start[EXONS$exon=="E7"] # before E7
intron7_start <- EXONS$end[EXONS$exon=="E7"]   # after E7
intron7_end   <- EXONS$start[EXONS$exon=="E8"] # before E8

message("Intron 6: chr5:", intron6_start, "-", intron6_end,
        " (", intron6_end-intron6_start, "bp)")
message("Intron 7: chr5:", intron7_start, "-", intron7_end,
        " (", intron7_end-intron7_start, "bp)")

# Load H3K27ac peaks
load_peaks <- function(file, label) {
  df <- read.table(gzfile(file), header=FALSE, sep='\t',
    col.names=c('chr','start','end','name','score',
                'strand','fc','pval','qval','summit'))
  df$label <- label
  df
}

peaks_ctrl1 <- load_peaks(
  'data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep1.narrowPeak.gz', 'CTRL_Rep1')
peaks_ctrl2 <- load_peaks(
  'data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep2.narrowPeak.gz', 'CTRL_Rep2')
peaks_vpa1  <- load_peaks(
  'data/external/h3k27ac_gse246399/H3K27ac_VPA_Rep1.narrowPeak.gz',  'VPA_Rep1')
all_peaks <- rbind(peaks_ctrl1, peaks_ctrl2, peaks_vpa1)

# Check intron 6
message('\n--- Intron 6 (E6-E7 boundary, flanks ASO target) ---')
int6 <- all_peaks[all_peaks$chr=='chr5' &
                  all_peaks$start >= intron6_start &
                  all_peaks$end   <= intron6_end, ]
if (nrow(int6) > 0) {
  message('PEAKS FOUND in intron 6:')
  print(int6[,c('chr','start','end','fc','qval','label')])
} else {
  message('No H3K27ac peaks in intron 6 in any condition')
}

# Check intron 7
message('\n--- Intron 7 (E7-E8 boundary, 3 prime end) ---')
int7 <- all_peaks[all_peaks$chr=='chr5' &
                  all_peaks$start >= intron7_start &
                  all_peaks$end   <= intron7_end, ]
if (nrow(int7) > 0) {
  message('PEAKS FOUND in intron 7:')
  print(int7[,c('chr','start','end','fc','qval','label')])
} else {
  message('No H3K27ac peaks in intron 7 in any condition')
}

# Extended check — 5kb around E7 (ASO target region)
message('\n--- 5kb window around E7 (ASO target ±5kb) ---')
e7_window <- all_peaks[all_peaks$chr=='chr5' &
                        all_peaks$start >= (EXONS$start[EXONS$exon=="E7"] - 5000) &
                        all_peaks$end   <= (EXONS$end[EXONS$exon=="E7"]   + 5000), ]
if (nrow(e7_window) > 0) {
  message('PEAKS near E7:')
  print(e7_window[,c('chr','start','end','fc','qval','label')])
} else {
  message('No H3K27ac peaks within 5kb of E7')
}

# Save intron results
intron_results <- data.frame(
  region = c('Intron 6', 'Intron 7', 'E7 ±5kb'),
  coordinates = c(
    paste0('chr5:', intron6_start, '-', intron6_end),
    paste0('chr5:', intron7_start, '-', intron7_end),
    paste0('chr5:', EXONS$start[EXONS$exon=="E7"]-5000, '-',
                    EXONS$end[EXONS$exon=="E7"]+5000)
  ),
  peaks_CTRL_Rep1 = c(
    sum(int6$label=='CTRL_Rep1'),
    sum(int7$label=='CTRL_Rep1'),
    sum(e7_window$label=='CTRL_Rep1')
  ),
  peaks_CTRL_Rep2 = c(
    sum(int6$label=='CTRL_Rep2'),
    sum(int7$label=='CTRL_Rep2'),
    sum(e7_window$label=='CTRL_Rep2')
  ),
  peaks_VPA_Rep1 = c(
    sum(int6$label=='VPA_Rep1'),
    sum(int7$label=='VPA_Rep1'),
    sum(e7_window$label=='VPA_Rep1')
  )
)
write.csv(intron_results,
  'results/smn2_enhancer/SMN2_intron67_H3K27ac_summary.csv',
  row.names=FALSE)
message('\nSaved: SMN2_intron67_H3K27ac_summary.csv')

# ── ITEM 10: Chr13 hotspot gene annotation ────────────────────────────────────
message('\n=== ITEM 10: Chromosome 13 DMR hotspot annotation ===')

# Chr13 hotspot region from meeting: ~60-80 Mb
CHR13_START <- 60000000
CHR13_END   <- 80000000

# Load ASO_VPA vs Scramble_CTRL DMRs (largest contrast)
dmr_file <- 'results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds'
if (file.exists(dmr_file)) {
  dmrs <- readRDS(dmr_file)
  message('Loaded DMRs: ', length(dmrs))

  # Filter to chr13 hotspot
  chr13_dmrs <- dmrs[seqnames(dmrs)=='chr13' &
                     start(dmrs) >= CHR13_START &
                     end(dmrs)   <= CHR13_END]
  message('Chr13 hotspot DMRs (60-80Mb): ', length(chr13_dmrs))

  if (length(chr13_dmrs) > 0) {
    # Annotate
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    anno <- annotatePeak(chr13_dmrs, tssRegion=c(-2000,2000),
      TxDb=txdb, annoDb='org.Hs.eg.db', verbose=FALSE)
    anno_df <- as.data.frame(anno)

    # Get unique genes
    genes <- unique(anno_df$SYMBOL[!is.na(anno_df$SYMBOL)])
    message('Unique genes in chr13 hotspot: ', length(genes))
    message('Top genes: ', paste(head(genes, 20), collapse=', '))

    write.csv(anno_df,
      'results/dmr_annotation/chr13_hotspot_annotated.csv',
      row.names=FALSE)

    # GO enrichment
    gene_ids <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
    if (length(gene_ids) >= 10) {
      go_res <- enrichGO(gene=gene_ids, OrgDb=org.Hs.eg.db,
        ont='BP', pAdjustMethod='BH',
        pvalueCutoff=0.05, readable=TRUE)

      if (!is.null(go_res) && nrow(go_res@result[go_res@result$p.adjust<0.05,]) > 0) {
        message('\nSignificant GO terms in chr13 hotspot:')
        print(head(go_res@result[go_res@result$p.adjust<0.05,
          c('Description','p.adjust','Count')], 10))
        write.csv(go_res@result,
          'results/dmr_annotation/chr13_hotspot_GO.csv',
          row.names=FALSE)

        pdf('results/dmr_annotation/chr13_hotspot_GO_dotplot.pdf',
            width=10, height=8)
        print(dotplot(go_res, showCategory=15,
          title='GO BP — Chr13 DMR hotspot (60-80Mb)\nASO_VPA vs Scramble_CTRL'))
        dev.off()
      } else {
        message('No significant GO terms in chr13 hotspot')
      }
    }

    # Also check chromatin/heterochromatin genes specifically
    message('\n--- Checking for heterochromatin/repeat genes ---')
    hetero_genes <- anno_df$SYMBOL[grep('KRAB|ZNF|LMNB|SATB|HP1|CBX|HDAC|DNMT|H3K|KDM|KMT',
      anno_df$SYMBOL, ignore.case=TRUE)]
    if (length(hetero_genes) > 0) {
      message('Chromatin-related genes: ', paste(unique(hetero_genes), collapse=', '))
    }

    # Region type summary
    message('\nAnnotation summary:')
    print(table(anno_df$annotation))
  }
} else {
  message('DMR file not found: ', dmr_file)
  # Try alternative
  dmr_files <- list.files('results/dmr', pattern='*.rds', full.names=TRUE)
  message('Available DMR files: ', paste(dmr_files, collapse=', '))
}

message('\nDone.')
