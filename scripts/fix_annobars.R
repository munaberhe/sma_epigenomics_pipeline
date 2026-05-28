.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomeInfoDb)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
OUT  <- 'results/dmr_annotation'

contrasts <- c('ASO_CTRL_vs_Scramble_CTRL',
               'ASO_VPA_vs_Scramble_CTRL',
               'Scramble_VPA_vs_Scramble_CTRL')

for (contrast in contrasts) {
  cat('Processing:', contrast, '\n')
  d  <- readRDS(paste0('results/dmr/dmr_', contrast, '.rds'))
  hc <- d[d$cytosinesCount >= 6]
  GenomeInfoDb::seqlevelsStyle(hc) <- 'UCSC'

  anno <- annotatePeak(hc, tssRegion=c(-3000,3000),
                       TxDb=txdb, annoDb='org.Hs.eg.db')

  cairo_pdf(file.path(OUT, paste0(contrast, '_annotation_bar.pdf')),
            width=10, height=6)
  print(plotAnnoBar(anno,
        title=paste0('Genomic Feature Distribution — ', contrast)))
  dev.off()
  cat(' done\n')
}
cat('Done.\n')
