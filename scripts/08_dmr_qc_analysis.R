# ============================================================
# DMRcaller genome-wide QC and supplementary analysis
# Muna Berhe · Zabet Lab · QMUL · 2026
# Covers:
#   1. Sample correlation / PCA
#   2. Permutation null distribution
#   3. DMR size distribution
#   4. CpG island overlap
#   5. Genome-wide methylation profile plots
# ============================================================

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

cat('Loading DMR objects...\n')
dmr_list <- lapply(CONTRASTS, function(ct) {
  d <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  d[d$cytosinesCount >= 6]
})
names(dmr_list) <- CONTRASTS

# ============================================================
# 1. Sample correlation heatmap and PCA
# ============================================================
cat('\n--- Analysis 1: Sample correlation / PCA ---\n')

prop_mat <- data.frame(
  ASO_CTRL      = dmr_list[['ASO_CTRL_vs_Scramble_CTRL']]$proportion1,
  Scramble_CTRL = dmr_list[['ASO_CTRL_vs_Scramble_CTRL']]$proportion2,
  ASO_VPA       = dmr_list[['ASO_VPA_vs_Scramble_CTRL']]$proportion1[
    seq_len(length(dmr_list[['ASO_CTRL_vs_Scramble_CTRL']]))],
  Scramble_VPA  = dmr_list[['Scramble_VPA_vs_Scramble_CTRL']]$proportion1[
    seq_len(length(dmr_list[['ASO_CTRL_vs_Scramble_CTRL']]))]
)
prop_mat <- prop_mat[complete.cases(prop_mat), ]

cor_mat <- cor(prop_mat, method='pearson')
cat('Pearson correlation matrix:\n')
print(round(cor_mat, 4))
write.csv(cor_mat, file.path(OUT, 'sample_correlation_matrix.csv'))

cor_melt <- melt(cor_mat)
pdf(file.path(OUT, 'sample_correlation_heatmap.pdf'), width=7, height=6)
print(ggplot(cor_melt, aes(Var1, Var2, fill=value)) +
  geom_tile(color='white') +
  geom_text(aes(label=round(value,3)), size=4, fontface='bold') +
  scale_fill_gradient2(low='#2166AC', mid='white', high='#B2182B',
                       midpoint=0.85, limits=c(0.7,1), name='Pearson r') +
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid=element_blank()) +
  labs(title='Sample methylation correlation',
       subtitle='Based on high-confidence DMRs (>=6 CpGs)', x='', y=''))
dev.off()

pca <- prcomp(t(prop_mat), scale.=TRUE)
pca_df <- data.frame(
  PC1       = pca$x[,1],
  PC2       = pca$x[,2],
  Condition = rownames(pca$x)
)
var_exp <- round(summary(pca)$importance[2,1:2]*100, 1)
cat('PC1 variance explained:', var_exp[1], '%\n')
cat('PC2 variance explained:', var_exp[2], '%\n')

pdf(file.path(OUT, 'sample_PCA.pdf'), width=7, height=6)
print(ggplot(pca_df, aes(PC1, PC2, label=Condition)) +
  geom_point(aes(color=Condition), size=5) +
  geom_text(vjust=-0.8, size=3.5) +
  scale_color_manual(values=c('ASO_CTRL'='#1B4F8A',
                               'Scramble_CTRL'='#6B7280',
                               'ASO_VPA'='#B2182B',
                               'Scramble_VPA'='#D97706')) +
  theme_minimal(base_size=13) +
  labs(title='PCA of condition methylation profiles',
       subtitle='Based on high-confidence DMRs (>=6 CpGs)',
       x=paste0('PC1 (', var_exp[1], '%)'),
       y=paste0('PC2 (', var_exp[2], '%)')) +
  theme(legend.position='none'))
dev.off()
cat('Correlation heatmap and PCA saved.\n')

# ============================================================
# 2. Permutation null distribution
# ============================================================
cat('\n--- Analysis 2: Permutation null distribution ---\n')

pdf(file.path(OUT, 'permutation_null_distribution.pdf'), width=12, height=5)
par(mfrow=c(1,3), mar=c(5,4,4,2))
for (i in seq_along(CONTRASTS)) {
  ct  <- CONTRASTS[i]
  d   <- dmr_list[[ct]]
  lab <- gsub('_vs_', ' vs\n', gsub('_', ' ', ct))
  obs_pvals <- d$pValue
  n <- min(length(obs_pvals), 5000)
  obs_sample <- sort(sample(obs_pvals, n))
  exp_uniform <- sort(runif(n))
  plot(-log10(exp_uniform), -log10(obs_sample),
       pch=20, cex=0.4,
       col=ifelse(obs_sample < 0.001, '#B2182B', '#6B7280'),
       main=lab, xlab='Expected -log10(p)', ylab='Observed -log10(p)',
       cex.main=0.9)
  abline(0, 1, col='black', lty=2, lwd=1.5)
  lambda <- median(qchisq(1-obs_pvals, 1)) / qchisq(0.5, 1)
  mtext(paste0('lambda = ', round(lambda, 3)), side=3, line=-1.5, cex=0.8)
}
dev.off()
cat('Permutation QQ plots saved.\n')

# ============================================================
# 3. DMR size distribution
# ============================================================
cat('\n--- Analysis 3: DMR size distribution ---\n')

size_df <- do.call(rbind, lapply(seq_along(CONTRASTS), function(i) {
  d <- dmr_list[[CONTRASTS[i]]]
  data.frame(
    width     = width(d),
    contrast  = gsub('_vs_Scramble_CTRL', '', CONTRASTS[i]),
    direction = ifelse(d$regionType == 'gain', 'Hypo', 'Hyper')
  )
}))

pdf(file.path(OUT, 'DMR_size_distribution.pdf'), width=10, height=5)
print(ggplot(size_df, aes(x=width, fill=direction)) +
  geom_histogram(bins=60, alpha=0.8, position='identity') +
  facet_wrap(~contrast, scales='free_y') +
  scale_fill_manual(values=c('Hypo'='#2166AC', 'Hyper'='#B2182B')) +
  scale_x_continuous(limits=c(0, 3000)) +
  theme_minimal(base_size=12) +
  labs(title='DMR size distribution by contrast',
       subtitle='High-confidence DMRs (>=6 CpGs); x-axis capped at 3kb',
       x='DMR width (bp)', y='Count', fill='Direction'))
dev.off()

cat('DMR width summary per contrast:\n')
for (ct in CONTRASTS) {
  d <- dmr_list[[ct]]
  cat(ct, '— median:', median(width(d)), 'bp | mean:', round(mean(width(d))),
      'bp | max:', max(width(d)), 'bp\n')
}
cat('DMR size distribution saved.\n')

# ============================================================
# 4. CpG island overlap
# ============================================================
cat('\n--- Analysis 4: CpG island overlap ---\n')

cpgi_file <- 'data/reference/cpg_islands_hg38.bed'
if (!file.exists(cpgi_file)) {
  dir.create(dirname(cpgi_file), showWarnings=FALSE, recursive=TRUE)
  cat('Downloading CpG islands from UCSC...\n')
  tryCatch({
    cpgi_url <- 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz'
    tmp <- tempfile()
    download.file(cpgi_url, tmp, quiet=TRUE)
    cpgi_raw <- read.table(gzfile(tmp), header=FALSE, sep='\t')
    cpgi_bed <- data.frame(
      chr   = cpgi_raw[,2],
      start = cpgi_raw[,3],
      end   = cpgi_raw[,4],
      name  = cpgi_raw[,5]
    )
    write.table(cpgi_bed, cpgi_file, sep='\t', quote=FALSE,
                row.names=FALSE, col.names=FALSE)
    cat('CpG islands downloaded and saved.\n')
  }, error=function(e) {
    cat('Download failed:', conditionMessage(e), '\n')
    cat('Skipping CpG island analysis.\n')
    cpgi_file <<- NULL
  })
}

if (!is.null(cpgi_file) && file.exists(cpgi_file)) {
  cpgi <- import(cpgi_file, format='BED')
  seqlevelsStyle(cpgi) <- 'UCSC'
  shores  <- c(flank(cpgi, 2000, start=TRUE), flank(cpgi, 2000, start=FALSE))
  shelves <- setdiff(c(flank(cpgi, 4000, start=TRUE),
                        flank(cpgi, 4000, start=FALSE)), shores)

  results <- do.call(rbind, lapply(seq_along(CONTRASTS), function(i) {
    d <- dmr_list[[CONTRASTS[i]]]
    seqlevelsStyle(d) <- 'UCSC'
    total    <- length(d)
    n_island <- length(subsetByOverlaps(d, cpgi))
    n_shore  <- length(subsetByOverlaps(d, shores))
    n_shelf  <- length(subsetByOverlaps(d, shelves))
    n_sea    <- total - n_island - n_shore - n_shelf
    data.frame(
      Contrast = gsub('_vs_Scramble_CTRL','',CONTRASTS[i]),
      Island   = round(n_island/total*100, 1),
      Shore    = round(n_shore/total*100, 1),
      Shelf    = round(n_shelf/total*100, 1),
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
    scale_fill_manual(values=c('Island'='#1B4F8A', 'Shore'='#4393C3',
                                'Shelf'='#92C5DE', 'Open_Sea'='#D1D5DB')) +
    theme_minimal(base_size=13) +
    labs(title='DMR distribution relative to CpG islands',
         x='', y='Percentage of DMRs (%)', fill='CpG context') +
    geom_text(aes(label=paste0(Percent,'%')),
              position=position_stack(vjust=0.5),
              size=3.5, color='white', fontface='bold'))
  dev.off()
  cat('CpG island overlap saved.\n')
}

# ============================================================
# 5. Genome-wide methylation distribution
# ============================================================
cat('\n--- Analysis 5: Genome-wide methylation distribution ---\n')

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
cat('Methylation distribution plot saved.\n')

# ============================================================
# Done
# ============================================================
cat('\n=== All analyses complete. Outputs in:', OUT, '===\n')
cat('Files generated:\n')
for (f in list.files(OUT)) cat(' -', f, '\n')
