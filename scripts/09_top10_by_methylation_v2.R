.libPaths('~/R/library')
suppressPackageStartupMessages({
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

  anno_df <- read.csv(file.path(OUT, paste0(contrast, '_annotated.csv')))

  # Filter to hypomethylated only
  hypo <- anno_df[anno_df$regionType == 'loss' & !is.na(anno_df$SYMBOL), ]

  # Composite filters:
  # 1. Minimum 10 CpGs per bin — high confidence estimate
  # 2. Baseline methylation >= 0.20 — meaningful biological drop
  # 3. Exclude chrX — likely karyotype artefact
  hypo <- hypo[hypo$cytosinesCount >= 10 &
               hypo$proportion2 >= 0.20 &
               hypo$seqnames != 'chrX', ]

  # Methylation difference (negative = hypo in treatment)
  hypo$meth_diff <- hypo$proportion1 - hypo$proportion2

  # Deduplicate — keep one row per gene, the one with largest meth_diff
  hypo <- hypo[order(-hypo$meth_diff), ]
  hypo <- hypo[!duplicated(hypo$SYMBOL), ]

  # Sort by meth_diff descending, p-value as tiebreaker
  hypo <- hypo[order(-hypo$meth_diff, hypo$pValue), ]

  top10 <- head(hypo[, c('seqnames','start','end','cytosinesCount',
                          'regionType','proportion1','proportion2',
                          'pValue','meth_diff','annotation',
                          'SYMBOL','GENENAME')], 10)

  outfile <- file.path(OUT, paste0(contrast, '_top10_hypo_v2.csv'))
  write.csv(top10, outfile, row.names=FALSE)

  cat('Top 10 hypomethylated (cytosinesCount>=10, proportion2>=0.20, no chrX):\n')
  print(top10[, c('SYMBOL','meth_diff','pValue','cytosinesCount','annotation')])
  cat('\n')
}
cat('Done.\n')
