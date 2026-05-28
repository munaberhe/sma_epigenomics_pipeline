.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(data.table)
})

OUT <- 'results/smn2_masked_profile'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

CHR5_DIR <- 'results/alignments_smn1_masked/chr5_cx'

SAMPLES <- list(
  ASO_CTRL      = c('ASO_CTRL_1', 'ASO_CTRL_2', 'ASO_CTRL_3'),
  ASO_VPA       = c('ASO_VPA_1',  'ASO_VPA_2',  'ASO_VPA_3'),
  Scramble_CTRL = c('Scramble_CTRL_1','Scramble_CTRL_2','Scramble_CTRL_3'),
  Scramble_VPA  = c('Scramble_VPA_1', 'Scramble_VPA_2', 'Scramble_VPA_3')
)

# SMN2 and SMN1 locus windows
SMN2_CHR   <- 'chr5'
SMN2_START <- 70040000
SMN2_END   <- 70090000
SMN1_START <- 70910000
SMN1_END   <- 70960000

cat('Loading and pooling chr5 CX reports...\n')

load_and_pool <- function(sample_names, condition_name) {
  cat('  Pooling', condition_name, '...\n')
  dts <- lapply(sample_names, function(s) {
    path <- file.path(CHR5_DIR, paste0(s, '_chr5.CX_report.txt'))
    dt <- fread(path, header=FALSE, sep='\t',
                col.names=c('chr','pos','strand','M','U','context','trinuc'),
                showProgress=FALSE)
    dt[context == 'CG']
  })
  # Pool by summing M and U
  merged <- Reduce(function(a, b) {
    m <- merge(a[, .(chr,pos,strand,M,U)],
               b[, .(chr,pos,strand,M,U)],
               by=c('chr','pos','strand'), suffixes=c('','_b'))
    m$M <- m$M + m$M_b
    m$U <- m$U + m$U_b
    m[, .(chr,pos,strand,M,U)]
  }, dts)
  merged$coverage <- merged$M + merged$U
  merged$methylation <- ifelse(merged$coverage > 0, merged$M/merged$coverage, NA)
  merged[coverage >= 3]
}

pooled <- lapply(names(SAMPLES), function(cond) {
  load_and_pool(SAMPLES[[cond]], cond)
})
names(pooled) <- names(SAMPLES)

# Convert to DMRcaller-compatible RLElist format using readBismark approach
# Build GRanges objects per condition for plotMethylation
make_gr <- function(dt, chr, start, end) {
  sub <- dt[chr == get('chr', envir=parent.frame()) &
            pos >= start & pos <= end]
  GRanges(
    seqnames = sub$chr,
    ranges   = IRanges(sub$pos, sub$pos),
    strand   = sub$strand,
    coverage = sub$coverage,
    methylation = sub$methylation
  )
}

# Use DMRcaller's plotMethylationData function
# Need to build the methylationData objects using readBismark

cat('Building methylation data objects...\n')

build_meth_data <- function(dt) {
  # DMRcaller expects a specific GRanges with readsM and readsN metadata
  gr <- GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(dt$pos, dt$pos),
    strand   = dt$strand
  )
  gr$readsM <- dt$M
  gr$readsN <- dt$coverage
  gr$methylation <- dt$methylation
  gr
}

meth_objects <- lapply(pooled, build_meth_data)

# Plot using DMRcaller plotMethylationData
contrasts_list <- list(
  list(g1='ASO_CTRL',    g2='Scramble_CTRL', label='ASO effect (CTRL background)'),
  list(g1='ASO_VPA',     g2='Scramble_VPA',  label='ASO effect (VPA background)'),
  list(g1='ASO_CTRL',    g2='ASO_VPA',       label='VPA effect (ASO background)')
)

for (ct in contrasts_list) {
  cat('Plotting', ct$label, '...\n')
  fname <- paste0('SMN2_masked_dmrcaller_', ct$g1, '_vs_', ct$g2, '.pdf')

  cairo_pdf(file.path(OUT, fname), width=12, height=8)
  par(mfrow=c(2,1), mar=c(4,4,3,2))

  for (gene in c('SMN1','SMN2')) {
    gs  <- if (gene=='SMN1') SMN1_START else SMN2_START
    ge  <- if (gene=='SMN1') SMN1_END   else SMN2_END

    d1 <- meth_objects[[ct$g1]]
    d2 <- meth_objects[[ct$g2]]

    # Subset to region
    d1_reg <- d1[seqnames(d1)==SMN2_CHR & start(d1)>=gs & start(d1)<=ge]
    d2_reg <- d2[seqnames(d2)==SMN2_CHR & start(d2)>=gs & start(d2)<=ge]

    if (length(d1_reg)==0 && length(d2_reg)==0) {
      plot.new(); title(paste(gene, 'no data'))
      next
    }

    # Merge positions
    pos1 <- start(d1_reg); m1 <- d1_reg$methylation
    pos2 <- start(d2_reg); m2 <- d2_reg$methylation

    all_pos <- sort(unique(c(pos1, pos2)))
    met1 <- m1[match(all_pos, pos1)]
    met2 <- m2[match(all_pos, pos2)]

    plot(all_pos, met1, type='p', pch=20, cex=0.8,
         col=ifelse(ct$g1 %in% c('ASO_CTRL','ASO_VPA'), '#1B4F8A', '#6B7280'),
         xlim=c(gs,ge), ylim=c(0,1.1),
         xlab=paste0('genomic coordinate on chromosome ', SMN2_CHR),
         ylab='methylation proportion',
         main=paste0(gene, ' ... ', ct$g1, ' vs ', ct$g2,
                     ' (', ct$label, ')'))
    points(all_pos, met2, pch=20, cex=0.8,
           col=ifelse(ct$g2 %in% c('ASO_CTRL','ASO_VPA'), '#1B4F8A', '#6B7280'))
    lines(all_pos[!is.na(met1)], met1[!is.na(met1)],
          col=ifelse(ct$g1 %in% c('ASO_CTRL','ASO_VPA'), '#1B4F8A', '#6B7280'),
          lwd=1.5)
    lines(all_pos[!is.na(met2)], met2[!is.na(met2)],
          col=ifelse(ct$g2 %in% c('Scramble_CTRL','Scramble_VPA'), '#D97706', '#B2182B'),
          lwd=1.5)
    legend('topright', legend=c(ct$g1, ct$g2),
           col=c('#1B4F8A','#D97706'), lwd=2, pch=20, bty='n')

    # Mark exon 7 for SMN2
    if (gene=='SMN2') {
      rect(70070697, 0, 70070817, 1.05, col=rgb(1,0,0,0.1), border='red')
      text((70070697+70070817)/2, 1.08, 'E7', col='red', cex=0.8, font=2)
    }
  }
  dev.off()
  cat(' Saved:', fname, '\n')
}

# Also save the weighted mean table
cat('\nWeighted mean methylation (masked, min 3x coverage):\n')
means <- do.call(rbind, lapply(names(pooled), function(cond) {
  dt <- pooled[[cond]]
  smn2 <- dt[chr=='chr5' & pos>=SMN2_START & pos<=SMN2_END & coverage>=3]
  smn1 <- dt[chr=='chr5' & pos>=SMN1_START & pos<=SMN1_END & coverage>=3]
  rbind(
    data.frame(condition=cond, locus='SMN2',
               wmean=round(weighted.mean(smn2$methylation, smn2$coverage, na.rm=TRUE),3),
               n_cpg=nrow(smn2)),
    data.frame(condition=cond, locus='SMN1',
               wmean=round(weighted.mean(smn1$methylation, smn1$coverage, na.rm=TRUE),3),
               n_cpg=nrow(smn1))
  )
}))
print(means[order(means$locus, means$condition),])
write.csv(means, file.path(OUT, 'SMN_masked_weighted_mean_methylation_v2.csv'),
          row.names=FALSE)

cat('\n=== Done ===\n')
for (f in list.files(OUT, pattern='\\.pdf$')) cat(' -', f, '\n')
