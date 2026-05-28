.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(reshape2)
  library(rtracklayer)
})

OUT <- 'results/dmr_qc'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

CONTRASTS <- c('ASO_CTRL_vs_Scramble_CTRL',
               'ASO_VPA_vs_Scramble_CTRL',
               'Scramble_VPA_vs_Scramble_CTRL')

dmr_list <- lapply(CONTRASTS, function(ct) {
  d <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  d[d$cytosinesCount >= 6]
})
names(dmr_list) <- CONTRASTS

# ============================================================
# 4. CpG island overlap — chunked for large contrasts
# ============================================================
cat('--- Analysis 4: CpG island overlap ---\n')

cpgi <- import('data/reference/cpg_islands_hg38.bed', format='BED')
seqlevelsStyle(cpgi) <- 'UCSC'
shores <- reduce(c(flank(cpgi, 2000, start=TRUE),
                    flank(cpgi, 2000, start=FALSE)))

overlap_pct <- function(dmrs, features) {
  # Process by chromosome to avoid memory issues
  chrs <- unique(as.character(seqnames(dmrs)))
  total <- 0
  for (chr in chrs) {
    d_chr <- dmrs[seqnames(dmrs)==chr]
    f_chr <- features[seqnames(features)==chr]
    if (length(f_chr) > 0)
      total <- total + length(subsetByOverlaps(d_chr, f_chr))
    else
      total <- total + 0
  }
  total
}

results <- do.call(rbind, lapply(seq_along(CONTRASTS), function(i) {
  d <- dmr_list[[CONTRASTS[i]]]
  seqlevelsStyle(d) <- 'UCSC'
  total    <- length(d)
  cat('Processing', CONTRASTS[i], '(', total, 'DMRs)...\n')
  n_island <- overlap_pct(d, cpgi)
  n_shore  <- overlap_pct(d, shores) - n_island
  n_shore  <- max(n_shore, 0)
  n_sea    <- total - n_island - n_shore
  cat('  Island:', n_island, '| Shore:', n_shore, '| Sea:', n_sea, '\n')
  data.frame(
    Contrast = gsub('_vs_Scramble_CTRL','',CONTRASTS[i]),
    Island   = round(n_island/total*100, 1),
    Shore    = round(n_shore/total*100, 1),
    Open_Sea = round(n_sea/total*100, 1)
  )
}))

cat('CpG context distribution (%):\n')
print(results)
write.csv(results, file.path(OUT, 'CpG_island_overlap.csv'), row.names=FALSE)

res_melt <- melt(results, id.vars='Contrast',
                 variable.name='Context', value.name='Percent')
pdf(file.path(OUT, 'CpG_island_overlap.pdf'), width=9, height=5)
print(ggplot(res_melt, aes(x=Contrast, y=Percent, fill=Context)) +
  geom_bar(stat='identity', position='stack', width=0.6) +
  scale_fill_manual(values=c('Island'='#1B4F8A',
                              'Shore'='#4393C3',
                              'Open_Sea'='#D1D5DB')) +
  theme_minimal(base_size=13) +
  labs(title='DMR distribution relative to CpG islands',
       x='', y='Percentage of DMRs (%)', fill='CpG context') +
  geom_text(aes(label=paste0(Percent,'%')),
            position=position_stack(vjust=0.5),
            size=3.5, color='white', fontface='bold'))
dev.off()
cat('CpG island overlap saved.\n')

# ============================================================
# 5. Genome-wide methylation distribution
# ============================================================
cat('--- Analysis 5: Methylation distribution ---\n')

meth_df <- do.call(rbind, lapply(seq_along(CONTRASTS), function(i) {
  d  <- dmr_list[[CONTRASTS[i]]]
  ct <- CONTRASTS[i]
  rbind(
    data.frame(methylation = d$proportion1,
               condition   = strsplit(ct,'_vs_')[[1]][1]),
    data.frame(methylation = d$proportion2,
               condition   = strsplit(ct,'_vs_')[[1]][2])
  )
}))
meth_df <- meth_df[!duplicated(paste(meth_df$condition,
                                      round(meth_df$methylation,3))),]

pdf(file.path(OUT, 'genome_wide_methylation_distribution.pdf'), width=9, height=5)
print(ggplot(meth_df, aes(x=methylation, color=condition, fill=condition)) +
  geom_density(alpha=0.15, size=0.9) +
  scale_color_manual(values=c('ASO_CTRL'='#1B4F8A',
                               'Scramble_CTRL'='#6B7280',
                               'ASO_VPA'='#B2182B',
                               'Scramble_VPA'='#D97706')) +
  scale_fill_manual(values=c('ASO_CTRL'='#1B4F8A',
                              'Scramble_CTRL'='#6B7280',
                              'ASO_VPA'='#B2182B',
                              'Scramble_VPA'='#D97706')) +
  theme_minimal(base_size=13) +
  labs(title='Genome-wide CpG methylation distribution',
       subtitle='Density across high-confidence DMR loci per condition',
       x='Methylation proportion', y='Density',
       color='Condition', fill='Condition') +
  xlim(0,1))
dev.off()
cat('Methylation distribution saved.\n')

cat('\n=== Done ===\n')
for (f in list.files(OUT)) cat(' -', f, '\n')
