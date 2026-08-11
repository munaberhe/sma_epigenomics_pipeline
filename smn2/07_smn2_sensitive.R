.libPaths(c('~/R/library', .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(dplyr)
})
setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/smn2_local_dmr'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

SMN2_START <- 70049638 - 10000
SMN2_END   <- 70078522 + 10000

CONTRASTS <- list(
  list(name='ASO_CTRL_vs_Scramble_CTRL',
       cond1=c('ASO_CTRL_1','ASO_CTRL_2','ASO_CTRL_3'),
       cond2=c('Scramble_CTRL_1','Scramble_CTRL_2','Scramble_CTRL_3'),
       masked=TRUE),
  list(name='Scramble_VPA_vs_Scramble_CTRL',
       cond1=c('Scramble_VPA_1','Scramble_VPA_2','Scramble_VPA_3'),
       cond2=c('Scramble_CTRL_1','Scramble_CTRL_2','Scramble_CTRL_3'),
       masked=TRUE),
  list(name='ASO_VPA_vs_ASO_CTRL',
       cond1=c('ASO_VPA_1','ASO_VPA_2','ASO_VPA_3'),
       cond2=c('ASO_CTRL_1','ASO_CTRL_2','ASO_CTRL_3'),
       masked=TRUE),
  list(name='ASO_VPA_vs_Scramble_VPA',
       cond1=c('ASO_VPA_1','ASO_VPA_2','ASO_VPA_3'),
       cond2=c('Scramble_VPA_1','Scramble_VPA_2','Scramble_VPA_3'),
       masked=TRUE)
)

read_masked <- function(samples) {
  grs <- lapply(samples, function(s) {
    f <- file.path('results/alignments/bs/by_chr',
                   paste0(s, '_chr5.CpG_report.txt.gz'))
    d <- read.table(gzfile(f), header=FALSE, sep='\t',
      col.names=c('chr','pos','strand','countM','countU','context','tri'),
      colClasses=c('character','integer','character','integer',
                   'integer','character','character'))
    d <- d[d$context=='CG' & d$pos >= SMN2_START & d$pos <= SMN2_END, ]
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos, d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

PARAM_SETS <- list(
  list(label="5pct_100bp", minDiff=0.05, binSize=100, minCyto=3, pval=0.05),
  list(label="2pct_50bp",  minDiff=0.02, binSize=50,  minCyto=2, pval=0.05),
  list(label="1pct_50bp",  minDiff=0.01, binSize=50,  minCyto=2, pval=0.10),
  # Radu: genome-wide 300bp window tuned on VPA signal may be too coarse for
  # the small ASO-only differences at SMN2 — test 200bp and 500bp here
  list(label="5pct_200bp", minDiff=0.05, binSize=200, minCyto=3, pval=0.05),
  list(label="5pct_500bp", minDiff=0.05, binSize=500, minCyto=3, pval=0.05)
)

all_results <- list()

for (ct in CONTRASTS) {
  message('\nContrast: ', ct$name)
  p1 <- read_masked(ct$cond1)
  p2 <- read_masked(ct$cond2)
  smn2_region <- GRanges('chr5', IRanges(SMN2_START, SMN2_END))

  for (ps in PARAM_SETS) {
    message('  params: ', ps$label)
    dmrs <- computeDMRs(p1, p2,
      regions=smn2_region, context='CG', method='bins',
      binSize=ps$binSize, minCytosinesCount=ps$minCyto,
      minReadsPerCytosine=3, pValueThreshold=ps$pval, test='score')
    if (length(dmrs) > 0) {
      md <- abs(dmrs$proportion1 - dmrs$proportion2)
      dmrs <- dmrs[md >= ps$minDiff]
    }
    message('    DMRs found: ', length(dmrs))
    if (length(dmrs) > 0) {
      df <- as.data.frame(dmrs)
      df$contrast  <- ct$name
      df$threshold <- ps$label
      df$near_E7   <- df$start >= 70070000 & df$end <= 70082000
      all_results[[paste(ct$name, ps$label, sep="_")]] <- df
      print(df[, c('seqnames','start','end','proportion1','proportion2',
                   'pValue','regionType','near_E7')])
    }
  }
}

non_empty <- Filter(function(x) nrow(x) > 0, all_results)
if (length(non_empty) > 0) {
  combined <- do.call(rbind, non_empty)
  write.csv(combined, file.path(OUT, 'SMN2_sensitive_DMRs_all_contrasts.csv'),
            row.names=FALSE)
  message('\nSaved: SMN2_sensitive_DMRs_all_contrasts.csv')
  message('Total DMRs found: ', nrow(combined))
  message('\nSummary by contrast and threshold:')
  print(table(combined$contrast, combined$threshold))
} else {
  message('No DMRs found at SMN2 in any contrast at any threshold')
  writeLines('No SMN2 DMRs found in any contrast',
             file.path(OUT, 'SMN2_sensitive_all_contrasts_result.txt'))
}
message('\nDone.')


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
