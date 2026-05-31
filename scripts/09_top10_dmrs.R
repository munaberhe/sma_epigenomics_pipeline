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
               'Scramble_VPA_vs_Scramble_CTRL',
               'ASO_VPA_vs_ASO_CTRL',
               'ASO_VPA_vs_Scramble_VPA')

for (contrast in contrasts) {
  cat('Processing:', contrast, '\n')

  anno_df <- read.csv(file.path(OUT, paste0(contrast, '_annotated.csv')))

  # regionType == 'gain': proportion1 < proportion2
  # meaning treatment is hypomethylated relative to reference
  hypo <- anno_df[anno_df$regionType == 'gain' & !is.na(anno_df$SYMBOL), ]

  # methylation difference — positive = bigger drop in treatment
  hypo$meth_diff <- hypo$proportion2 - hypo$proportion1

  # Stage 1: filter to meaningful methylation difference
  # 0.20 is already the DMRcaller minDiff threshold but making it
  # explicit here ensures any annotation artefacts are excluded
  hypo <- hypo[hypo$meth_diff >= 0.20, ]

  # Stage 2: deduplicate by gene, keeping the DMR with largest effect
  hypo <- hypo[order(-hypo$meth_diff, hypo$pValue), ]
  hypo <- hypo[!duplicated(hypo$SYMBOL), ]

  # Stage 3: rank by p-value — most statistically confident first
  hypo <- hypo[order(hypo$pValue), ]

  top10 <- head(hypo[, c('seqnames','start','end','cytosinesCount',
                          'regionType','proportion1','proportion2',
                          'pValue','meth_diff','annotation',
                          'SYMBOL','GENENAME')], 10)

  outfile <- file.path(OUT, paste0(contrast, '_top10_hypo_v2.csv'))
  write.csv(top10, outfile, row.names=FALSE)

  cat('Top 10 hypomethylated (meth_diff>=0.20, ranked by p-value):\n')
  print(top10[, c('SYMBOL','meth_diff','pValue','cytosinesCount','annotation')])
  cat('\n')
}
cat('Done.\n')
