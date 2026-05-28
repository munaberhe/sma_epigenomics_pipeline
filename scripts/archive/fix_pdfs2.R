.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
OUT  <- 'results/dmr_annotation'

contrasts <- c('ASO_VPA_vs_Scramble_CTRL',
               'Scramble_VPA_vs_Scramble_CTRL')

for (contrast in contrasts) {
  cat('Regenerating annotation bar for:', contrast, '\n')
  anno_df <- read.csv(file.path(OUT, paste0(contrast, '_annotated.csv')))

  gr <- GRanges(
    seqnames = anno_df$seqnames,
    ranges   = IRanges(anno_df$start, anno_df$end)
  )
  mcols(gr) <- anno_df[, !names(anno_df) %in% c('seqnames','start','end')]

  anno <- annotatePeak(gr, tssRegion=c(-3000,3000),
                       TxDb=txdb, annoDb='org.Hs.eg.db')

  cairo_pdf(file.path(OUT, paste0(contrast, '_annotation_bar.pdf')),
            width=10, height=6)
  print(plotAnnoBar(anno,
        title=paste0('Genomic Feature Distribution — ', contrast)))
  dev.off()
  cat(' done\n')
}
cat('Done.\n')
