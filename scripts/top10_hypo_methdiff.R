.libPaths('~/R/library')
library(dplyr)
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

contrasts <- c('ASO_CTRL_vs_Scramble_CTRL',
               'Scramble_VPA_vs_Scramble_CTRL',
               'ASO_VPA_vs_Scramble_CTRL',
               'ASO_VPA_vs_ASO_CTRL',
               'ASO_VPA_vs_Scramble_VPA')

all_hypo <- lapply(contrasts, function(ct) {
  f <- paste0('results/dmr_annotation/', ct, '_annotated.csv')
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, stringsAsFactors=FALSE)
  df %>%
    filter(regionType=='gain', pValue < 0.01,
           sumReadsN1 >= 10, sumReadsN2 >= 10) %>%
    mutate(meth_diff = proportion2 - proportion1,
           contrast = ct) %>%
    filter(meth_diff >= 0.20) %>%
    arrange(desc(meth_diff)) %>%
    select(seqnames, start, end, proportion1, proportion2,
           meth_diff, pValue, annotation, SYMBOL, contrast) %>%
    head(5)
}) 
all_hypo <- bind_rows(all_hypo)
all_hypo <- all_hypo %>% arrange(desc(meth_diff))
print(all_hypo, row.names=FALSE)
write.csv(all_hypo, 'results/dmr_annotation/top_hypo_by_methdiff.csv', row.names=FALSE)
message("Saved: top_hypo_by_methdiff.csv")
