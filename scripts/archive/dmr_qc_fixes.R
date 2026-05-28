.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ggplot2)
  library(reshape2)
  library(data.table)
})

OUT <- 'results/dmr_qc'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ============================================================
# FIX 1 — PCA with corrected labels
# ============================================================
cat('--- Fix 1: PCA with corrected labels ---\n')

CONTRASTS <- c('ASO_CTRL_vs_Scramble_CTRL',
               'ASO_VPA_vs_Scramble_CTRL',
               'Scramble_VPA_vs_Scramble_CTRL')

dmr_list <- lapply(CONTRASTS, function(ct) {
  d <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  d[d$cytosinesCount >= 6]
})
names(dmr_list) <- CONTRASTS

prop_mat <- data.frame(
  ASO_CTRL      = dmr_list[['ASO_CTRL_vs_Scramble_CTRL']]$proportion1,
  Scramble_CTRL = dmr_list[['ASO_CTRL_vs_Scramble_CTRL']]$proportion2,
  ASO_VPA       = dmr_list[['ASO_VPA_vs_Scramble_CTRL']]$proportion1[
    seq_len(length(dmr_list[['ASO_CTRL_vs_Scramble_CTRL']]))],
  Scramble_VPA  = dmr_list[['Scramble_VPA_vs_Scramble_CTRL']]$proportion1[
    seq_len(length(dmr_list[['ASO_CTRL_vs_Scramble_CTRL']]))]
)
prop_mat <- prop_mat[complete.cases(prop_mat), ]

pca     <- prcomp(t(prop_mat), scale.=TRUE)
var_exp <- round(summary(pca)$importance[2,1:2]*100, 1)

pca_df <- data.frame(
  PC1       = pca$x[,1],
  PC2       = pca$x[,2],
  Condition = rownames(pca$x),
  VPA       = grepl('VPA', rownames(pca$x))
)

pdf(file.path(OUT, 'sample_PCA.pdf'), width=7, height=6)
print(ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(color=Condition), size=6) +
  geom_label(aes(label=Condition, color=Condition),
             vjust=-0.9, size=3.5, fontface='bold',
             fill='white', label.size=0.3, show.legend=FALSE) +
  scale_color_manual(values=c('ASO_CTRL'='#1B4F8A',
                               'Scramble_CTRL'='#6B7280',
                               'ASO_VPA'='#B2182B',
                               'Scramble_VPA'='#D97706')) +
  theme_minimal(base_size=13) +
  theme(legend.position='none') +
  expand_limits(y=c(min(pca_df$PC2)-8, max(pca_df$PC2)+12)) +
  labs(title='PCA of condition methylation profiles',
       subtitle=paste0('PC1 separates VPA vs non-VPA; PC2 separates ASO vs Scramble\n',
                       'Based on high-confidence DMRs (>=6 CpGs)'),
       x=paste0('PC1 (', var_exp[1], '%)'),
       y=paste0('PC2 (', var_exp[2], '%)')))
dev.off()
cat('PCA fixed.\n')

# ============================================================
# FIX 2 — Correlation heatmap from CX reports (chr1 only for speed)
# ============================================================
cat('--- Fix 2: Correlation heatmap from CX reports ---\n')

# Use chr1 from completed samples for genome-wide correlation
# Completed: ASO_CTRL_1, Scramble_CTRL_1/2/3, Scramble_VPA_1
cx_dir  <- 'results/alignments_smn1_masked/cx_report'
cx_files <- list(
  ASO_CTRL_1    = file.path(cx_dir, 'ASO_CTRL_1.CX_report.txt.gz'),
  Scramble_CTRL_1 = file.path(cx_dir, 'Scramble_CTRL_1.CX_report.txt.gz'),
  Scramble_CTRL_2 = file.path(cx_dir, 'Scramble_CTRL_2.CX_report.txt.gz'),
  Scramble_CTRL_3 = file.path(cx_dir, 'Scramble_CTRL_3.CX_report.txt.gz'),
  Scramble_VPA_1  = file.path(cx_dir, 'Scramble_VPA_1.CX_report.txt.gz')
)

cat('Loading CX reports (chr1 CpG only)...\n')
load_chr1_cpg <- function(path, name) {
  cat(' Reading', name, '...\n')
  dt <- fread(cmd=paste0('zcat ', path,
    ' | awk \'$1=="1" && $4+$5>=5 && $6=="CG"\''),
    header=FALSE, sep='\t',
    col.names=c('chr','pos','strand','M','U','context','trinuc'))
  dt$meth <- dt$M / (dt$M + dt$U)
  dt$key  <- paste0(dt$chr,'_',dt$pos,'_',dt$strand)
  dt[, c('key','meth'), with=FALSE]
}

meth_list <- mapply(load_chr1_cpg, cx_files, names(cx_files),
                    SIMPLIFY=FALSE)

# Merge on common CpGs
cat('Merging samples...\n')
merged <- Reduce(function(a,b) merge(a, b, by='key'), meth_list)
colnames(merged) <- c('key', names(cx_files))
merged <- merged[complete.cases(merged),]
cat('Common CpGs on chr1 (>=5x):', nrow(merged), '\n')

# Subsample for speed if very large
if (nrow(merged) > 100000) {
  set.seed(42)
  merged <- merged[sample(nrow(merged), 100000),]
  cat('Subsampled to 100,000 CpGs\n')
}

cor_mat <- cor(merged[,-1], method='pearson')
cat('Pearson correlation matrix:\n')
print(round(cor_mat, 4))
write.csv(cor_mat,
          file.path(OUT, 'sample_correlation_matrix_CX.csv'))

cor_melt <- melt(cor_mat)
pdf(file.path(OUT, 'sample_correlation_heatmap.pdf'), width=8, height=7)
print(ggplot(cor_melt, aes(Var1, Var2, fill=value)) +
  geom_tile(color='white', size=0.5) +
  geom_text(aes(label=round(value,3)), size=3.8, fontface='bold') +
  scale_fill_gradient2(low='#2166AC', mid='white', high='#B2182B',
                       midpoint=0.9, limits=c(0.7,1),
                       name='Pearson r') +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid=element_blank()) +
  labs(title='Sample CpG methylation correlation',
       subtitle='chr1 CpGs with >=5x coverage (5 of 12 samples with completed CX reports)',
       x='', y=''))
dev.off()
cat('Correlation heatmap fixed.\n')

# ============================================================
# FIX 3 — Methylation distribution from CX reports
# ============================================================
cat('--- Fix 3: Methylation distribution from CX reports ---\n')

# Use the merged chr1 data — assign condition labels
# Map sample to condition
sample_conditions <- c(
  ASO_CTRL_1      = 'ASO_CTRL',
  Scramble_CTRL_1 = 'Scramble_CTRL',
  Scramble_CTRL_2 = 'Scramble_CTRL',
  Scramble_CTRL_3 = 'Scramble_CTRL',
  Scramble_VPA_1  = 'Scramble_VPA'
)

meth_long <- do.call(rbind, lapply(names(sample_conditions), function(s) {
  data.frame(
    methylation = merged[[s]],
    sample      = s,
    condition   = sample_conditions[s]
  )
}))

condition_cols <- c('ASO_CTRL'='#1B4F8A',
                    'Scramble_CTRL'='#6B7280',
                    'Scramble_VPA'='#D97706')

pdf(file.path(OUT, 'genome_wide_methylation_distribution.pdf'), width=9, height=5)
print(ggplot(meth_long, aes(x=methylation, color=condition,
                             group=sample, linetype=sample)) +
  geom_density(alpha=0, size=0.7) +
  scale_color_manual(values=condition_cols) +
  scale_linetype_manual(values=rep(c('solid','dashed','dotted'), 5)) +
  theme_minimal(base_size=13) +
  guides(linetype='none') +
  labs(title='CpG methylation distribution per sample',
       subtitle='chr1 CpGs >=5x coverage | 5 samples with completed CX reports\nASO_VPA not yet available (extraction in progress)',
       x='Methylation proportion', y='Density', color='Condition') +
  xlim(0,1))
dev.off()
cat('Methylation distribution fixed.\n')

cat('\n=== All fixes done ===\n')
for (f in list.files(OUT)) cat(' -', f, '\n')
