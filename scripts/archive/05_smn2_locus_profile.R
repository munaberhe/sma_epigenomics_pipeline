.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(data.table)
})

OUT <- 'results/smn2_masked_profile'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

CHR5_DIR <- 'results/alignments_smn1_masked/chr5_cx'

# SMN locus coordinates on hg38
SMN2_START <- 70040000
SMN2_END   <- 70090000
SMN1_START <- 70910000
SMN1_END   <- 70960000

SAMPLES <- list(
  ASO_CTRL    = c('ASO_CTRL_1', 'ASO_CTRL_2', 'ASO_CTRL_3'),
  ASO_VPA     = c('ASO_VPA_1',  'ASO_VPA_2',  'ASO_VPA_3'),
  Scramble_CTRL = c('Scramble_CTRL_1', 'Scramble_CTRL_2', 'Scramble_CTRL_3'),
  Scramble_VPA  = c('Scramble_VPA_1',  'Scramble_VPA_2',  'Scramble_VPA_3')
)

COLOURS <- c(
  ASO_CTRL      = '#1B4F8A',
  Scramble_CTRL = '#6B7280',
  ASO_VPA       = '#B2182B',
  Scramble_VPA  = '#D97706'
)

cat('Loading chr5 CX reports...\n')

load_cx <- function(sample) {
  path <- file.path(CHR5_DIR, paste0(sample, '_chr5.CX_report.txt'))
  cat(' Reading', sample, '...\n')
  dt <- fread(path, header=FALSE, sep='\t',
              col.names=c('chr','pos','strand','M','U','context','trinuc'),
              showProgress=FALSE)
  dt <- dt[context == 'CG' & (M + U) >= 1]
  dt$meth <- dt$M / (dt$M + dt$U)
  dt
}

# Load all samples
all_data <- lapply(unlist(SAMPLES), load_cx)
names(all_data) <- unlist(SAMPLES)

# Pool replicates per condition
cat('Pooling replicates...\n')
pool_condition <- function(sample_names) {
  dts <- lapply(sample_names, function(s) all_data[[s]][, .(chr, pos, M, U)])
  merged <- Reduce(function(a, b) {
    merge(a, b, by=c('chr','pos'), suffixes=c('','_b'))
  }, dts)
  # Sum M and U across replicates
  m_cols <- grep('^M', names(merged), value=TRUE)
  u_cols <- grep('^U', names(merged), value=TRUE)
  merged$M_pool <- rowSums(merged[, ..m_cols])
  merged$U_pool <- rowSums(merged[, ..u_cols])
  merged$meth   <- merged$M_pool / (merged$M_pool + merged$U_pool)
  merged$cov    <- merged$M_pool + merged$U_pool
  merged[cov >= 3, .(chr, pos, meth, cov)]
}

pooled <- lapply(names(SAMPLES), function(cond) {
  cat('Pooling', cond, '...\n')
  dt <- pool_condition(SAMPLES[[cond]])
  dt$condition <- cond
  dt
})
pooled_all <- rbindlist(pooled)

# ── Plot function ─────────────────────────────────────────────────────────────
plot_locus <- function(data, start, end, gene_label, title, outfile) {
  d <- data[pos >= start & pos <= end]
  if (nrow(d) == 0) { cat('No data in region', gene_label, '\n'); return() }

  # Smooth with rolling mean per condition
  d <- d[order(condition, pos)]
  d[, meth_smooth := frollmean(meth, n=5, fill=NA, align='center'), by=condition]

  # Exon 7 coordinates for SMN2
  ex7_start <- 70070697
  ex7_end   <- 70070817

  p <- ggplot(d, aes(x=pos, y=meth, colour=condition)) +
    geom_point(alpha=0.3, size=0.8) +
    geom_line(aes(y=meth_smooth), size=0.9, na.rm=TRUE) +
    scale_colour_manual(values=COLOURS) +
    scale_y_continuous(limits=c(0,1), labels=scales::percent) +
    scale_x_continuous(labels=function(x) paste0(round(x/1e6,3), 'Mb')) +
    theme_minimal(base_size=12) +
    theme(legend.position='bottom',
          panel.grid.minor=element_blank()) +
    labs(title=title,
         subtitle=paste0('Masked alignment — chr5:', start, '-', end),
         x='Genomic coordinate (chr5)',
         y='CpG methylation proportion',
         colour='Condition')

  # Add exon 7 annotation for SMN2 region
  if (gene_label == 'SMN2') {
    p <- p +
      annotate('rect', xmin=ex7_start, xmax=ex7_end,
               ymin=-Inf, ymax=Inf,
               fill='red', alpha=0.1) +
      annotate('text', x=(ex7_start+ex7_end)/2, y=0.05,
               label='E7', colour='red', size=3.5, fontface='bold')
  }

  cairo_pdf(outfile, width=12, height=6)
  print(p)
  dev.off()
  cat('Saved:', outfile, '\n')
}

# ── SMN2 locus ───────────────────────────────────────────────────────────────
cat('\nPlotting SMN2 locus...\n')
plot_locus(pooled_all, SMN2_START, SMN2_END, 'SMN2',
           'SMN2 locus CpG methylation — masked alignment (all conditions)',
           file.path(OUT, 'SMN2_masked_all_conditions.pdf'))

# Individual contrasts
contrasts <- list(
  list(c1='ASO_CTRL', c2='Scramble_CTRL', label='ASO effect (CTRL background)'),
  list(c1='ASO_VPA',  c2='Scramble_VPA',  label='ASO effect (VPA background)'),
  list(c1='ASO_CTRL', c2='ASO_VPA',       label='VPA effect (ASO background)')
)

for (ct in contrasts) {
  d_sub <- pooled_all[condition %in% c(ct$c1, ct$c2)]
  plot_locus(d_sub, SMN2_START, SMN2_END, 'SMN2',
             paste0('SMN2 masked — ', ct$label),
             file.path(OUT, paste0('SMN2_masked_', ct$c1, '_vs_', ct$c2, '.pdf')))
}

# ── Weighted mean methylation table ─────────────────────────────────────────
cat('\nComputing weighted mean methylation at SMN loci...\n')

smn_means <- pooled_all[pos >= SMN2_START & pos <= SMN2_END,
  .(wmean = weighted.mean(meth, cov), n_cpg = .N),
  by=condition]
smn_means$locus <- 'SMN2'

smn1_means <- pooled_all[pos >= SMN1_START & pos <= SMN1_END,
  .(wmean = weighted.mean(meth, cov), n_cpg = .N),
  by=condition]
smn1_means$locus <- 'SMN1'

means_table <- rbind(smn_means, smn1_means)
means_table$wmean <- round(means_table$wmean, 3)
write.csv(means_table,
          file.path(OUT, 'SMN_masked_weighted_mean_methylation.csv'),
          row.names=FALSE)
cat('\nWeighted mean methylation (masked):\n')
print(means_table[order(locus, condition)])

cat('\n=== Done. Outputs in:', OUT, '===\n')
for (f in list.files(OUT)) cat(' -', f, '\n')
