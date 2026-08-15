.libPaths('~/R/library')
library(dplyr)
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

df <- read.csv('results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv')

hypo <- df %>%
  filter(regionType=='gain', pValue < 0.01,
         sumReadsN1 >= 10, sumReadsN2 >= 10) %>%
  mutate(meth_diff = proportion2 - proportion1) %>%
  filter(meth_diff >= 0.20) %>%
  arrange(desc(meth_diff)) %>%
  select(seqnames, start, end, proportion1, proportion2,
         meth_diff, pValue, annotation, SYMBOL, GENENAME) %>%
  head(20)

write.csv(hypo, 'results/dmr_annotation/top20_hypo_ASO_by_methdiff.csv', row.names=FALSE)
message("Saved: top20_hypo_ASO_by_methdiff.csv")
print(as.data.frame(hypo), row.names=FALSE)

df <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv')

hypo <- df %>%
  filter(regionType=='gain', pValue < 0.05,
         sumReadsN1 >= 10, sumReadsN2 >= 10) %>%
  mutate(meth_diff = proportion2 - proportion1) %>%
  filter(meth_diff >= 0.20) %>%
  arrange(desc(meth_diff)) %>%
  select(seqnames, start, end, proportion1, proportion2,
         meth_diff, pValue, annotation, SYMBOL, GENENAME) %>%
  head(20)

write.csv(hypo,
  'results/dmr_annotation/top20_hypo_ASO_VPA_vs_Scramble_CTRL.csv',
  row.names=FALSE)
message("Saved top20_hypo_ASO_VPA_vs_Scramble_CTRL.csv")
print(as.data.frame(hypo), row.names=FALSE)
