.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomeInfoDb)
  library(DMRcaller)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
OUT  <- 'results/dmr_annotation'

contrasts <- c('ASO_CTRL_vs_Scramble_CTRL',
               'ASO_VPA_vs_Scramble_CTRL',
               'Scramble_VPA_vs_Scramble_CTRL')

for (contrast in contrasts) {
  cat('Processing:', contrast, '\n')

  anno_df <- read.csv(file.path(OUT, paste0(contrast, '_annotated.csv')))

  # Filter to hypomethylated only (gain = proportion1 > proportion2
  # meaning treatment is more methylated... wait, need to check direction)
  # regionType == 'gain' means group1 proportion > group2 proportion
  # For VPA contrasts: group1=VPA, group2=Scramble_CTRL
  # gain = VPA MORE methylated than Scramble = hypermethylation in VPA
  # loss = VPA LESS methylated than Scramble = hypomethylation in VPA
  # So for hypomethylated in treatment: use regionType == 'loss'

  hypo <- anno_df[anno_df$regionType == "loss" & anno_df$proportion2 >= 0.05 & !is.na(anno_df$SYMBOL), ]

  # Methylation difference = proportion2 - proportion1 (positive = hypo in treatment)
  hypo$meth_diff <- hypo$proportion2 - hypo$proportion1

  # Rank by methylation difference (largest hypo change first)
  hypo <- hypo[order(-hypo$meth_diff), ]

  top10 <- head(hypo[, c('seqnames','start','end','cytosinesCount',
                          'regionType','proportion1','proportion2',
                          'pValue','meth_diff','annotation',
                          'SYMBOL','GENENAME')], 10)

  outfile <- file.path(OUT, paste0(contrast, '_top10_hypo_by_methdiff.csv'))
  write.csv(top10, outfile, row.names=FALSE)

  cat('Top 10 hypomethylated by methylation difference:\n')
  print(top10[, c('SYMBOL','meth_diff','pValue','annotation')])
  cat('\n')
}
