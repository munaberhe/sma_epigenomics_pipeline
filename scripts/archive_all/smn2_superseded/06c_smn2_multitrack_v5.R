#!/usr/bin/env Rscript
# 06c_smn2_multitrack_v5.R
# Multi-track SMN2 3-prime enhancer locus plot:
#   Track 1 -- H3K9me2 bigWig (CTRL vs ASO, kpPlotBigWig smooth signal)
#   Track 2 -- H3K27ac peaks (CTRL vs VPA, kpPlotRegions rectangles)
#   Track 3 -- Methylation lines (4 conditions, kpLines + kpPoints)
#   Track 4 -- CpG islands / H9 enhancers / ENCODE cCREs
#   Track 5 -- Gene model (SMN2 exon structure)
# Muna Berhe -- bt25018 -- QMUL MSc Bioinformatics

.libPaths('/data/home/bt25018/R/library')
suppressPackageStartupMessages({
  library(karyoploteR)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(annotatr)
  library(data.table)
})
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

OUT_DIR <- 'results/smn2_enhancer'
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# Locked palette
PAL <- c(ASO_CTRL='#1F3A5F', Scramble_CTRL='#6B7280',
         ASO_VPA='#C0392B', Scramble_VPA='#D4A017')
H3K9ME2_CTRL_COL <- '#4E9EC7'
H3K9ME2_ASO_COL  <- '#E67E22'
H3K27AC_CTRL_COL <- '#E31A1C'
H3K27AC_VPA_COL  <- '#8E44AD'

# Window -- DMR centred, SMN2 3-prime visible
ZOOM_CHR   <- 'chr5'
ZOOM_START <- 70044000
ZOOM_END   <- 70103000
zoom    <- toGRanges(data.frame(chr=ZOOM_CHR, start=ZOOM_START, end=ZOOM_END))
zoom_gr <- GRanges(ZOOM_CHR, IRanges(ZOOM_START, ZOOM_END))
DMR_POS <- 70088223

message(sprintf('Window: %s:%d-%d (%.1f kb)', ZOOM_CHR, ZOOM_START, ZOOM_END,
                (ZOOM_END-ZOOM_START)/1000))

# External data paths
H3K9ME2_BW <- list(
  CTRL_R1 = 'data/external/GSE167762_H3K9me2/GSM6063702_CTRvsInp_CTR_R1.bw',
  CTRL_R2 = 'data/external/GSE167762_H3K9me2/GSM6063706_CTRvsInp_CTR_R2.bw',
  ASO_R1  = 'data/external/GSE167762_H3K9me2/GSM6063703_ASOvsInp_ASO_R1.bw',
  ASO_R2  = 'data/external/GSE167762_H3K9me2/GSM6063707_ASOvsInp_ASO_R2.bw'
)
H3K27AC_NP <- list(
  CTRL = c('data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep1.narrowPeak.gz',
           'data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep2.narrowPeak.gz'),
  VPA  = c('data/external/h3k27ac_gse246399/H3K27ac_VPA_Rep1.narrowPeak.gz')
)

# Load H3K27ac peaks in window
load_narrowpeak <- function(files) {
  rows <- lapply(files[file.exists(files)], function(f) {
    d <- read.table(gzfile(f), header=FALSE, sep='\t',
                    col.names=c('chr','start','end','name','score','strand',
                                'fc','neglog10p','neglog10q','summit'))
    d[d$chr==ZOOM_CHR & d$end>=ZOOM_START & d$start<=ZOOM_END, ]
  })
  df <- do.call(rbind, rows)
  if (is.null(df) || nrow(df)==0) return(GRanges())
  # merge overlapping peaks
  gr <- reduce(GRanges(df$chr, IRanges(df$start, df$end)))
  gr
}
message('Loading H3K27ac peaks...')
h3k27ac_ctrl <- load_narrowpeak(H3K27AC_NP$CTRL)
h3k27ac_vpa  <- load_narrowpeak(H3K27AC_NP$VPA)
message('  CTRL peaks: ', length(h3k27ac_ctrl), '  VPA peaks: ', length(h3k27ac_vpa))

# Load annotation tracks
message('Loading CpG islands...')
cpg_anns <- build_annotations(genome='hg38', annotations='hg38_cpg_islands')
cgi <- subsetByOverlaps(cpg_anns, zoom_gr)
message('  CpG islands: ', length(cgi))

message('Loading H9 enhancers...')
enh_df <- read.table(gzfile('data/reference/H9_predicted_non_promoter_non_fragments.bed.gz'),
                     header=TRUE, sep='\t', stringsAsFactors=FALSE)
enh <- subsetByOverlaps(GRanges(enh_df$seqnames, IRanges(enh_df$start, enh_df$end)), zoom_gr)
message('  H9 enhancers: ', length(enh))

message('Loading ENCODE cCREs...')
ccre_df <- read.table('data/reference/encode_cCREs_hg38.bed', header=FALSE, sep='\t',
                       col.names=c('chr','start','end','id1','id2','type'))
ccre_sub <- ccre_df[grepl('ELS', ccre_df$type) & ccre_df$chr==ZOOM_CHR &
                    ccre_df$end>=ZOOM_START & ccre_df$start<=ZOOM_END, ]
ccre_gr <- if (nrow(ccre_sub)>0)
  GRanges(ccre_sub$chr, IRanges(ccre_sub$start, ccre_sub$end)) else GRanges()
message('  ENCODE cCREs: ', length(ccre_gr))

# Load methylation from masked CX files
CONDITIONS <- c('ASO_CTRL','Scramble_CTRL','ASO_VPA','Scramble_VPA')
BY_CHR <- 'results/alignments_smn1_masked/chr5_cx'

load_meth <- function(condition) {
  files <- list.files(BY_CHR,
    pattern=sprintf('^%s_[0-9]+_chr5\\.CX_report\\.txt$', condition),
    full.names=TRUE)
  if (length(files)==0) return(NULL)
  message('  loading: ', condition, ' (', length(files), ' reps)')
  df_list <- lapply(files, function(f) {
    d <- read.table(f, header=FALSE, sep='\t',
      col.names=c('chr','pos','strand','mC','uC','ctx','tri'))
    d <- d[d$ctx=='CG' & d$pos>=ZOOM_START & d$pos<=ZOOM_END, ]
    d$beta <- d$mC/(d$mC+d$uC)
    d[d$mC+d$uC>=3, c('pos','beta')]
  })
  merged <- do.call(rbind, df_list)
  if (is.null(merged) || nrow(merged)==0) return(NULL)
  merged <- aggregate(beta~pos, merged, mean, na.rm=TRUE)
  GRanges(ZOOM_CHR, IRanges(merged$pos, width=1), beta=merged$beta[order(merged$pos)])
}

message('Loading methylation profiles...')
meth <- setNames(lapply(CONDITIONS, load_meth), CONDITIONS)

rollmean_vec <- function(x, k=21) stats::filter(x, rep(1/k,k), sides=2)

# Build gene model
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

draw_genes <- function(kp, r0, r1) {
  all_genes <- genes(txdb, single.strand.genes.only=FALSE)
  if (is(all_genes, 'GRangesList')) all_genes <- unlist(all_genes)
  win_genes <- subsetByOverlaps(all_genes, zoom_gr)
  if (length(win_genes)==0) return(invisible(NULL))
  gene_ids <- names(win_genes)
  gene_names <- tryCatch(
    mapIds(org.Hs.eg.db, keys=gene_ids, column='SYMBOL',
           keytype='ENTREZID', multiVals='first'),
    error=function(e) setNames(gene_ids, gene_ids))
  gene_names[is.na(gene_names)] <- gene_ids[is.na(gene_names)]
  kpSegments(kp, chr=as.character(seqnames(win_genes)),
             x0=start(win_genes), x1=end(win_genes),
             y0=0.5, y1=0.5, col='grey30', lwd=1.5, r0=r0, r1=r1)
  ex_by_gene <- exonsBy(txdb, by='gene')
  ex_in_win  <- ex_by_gene[names(ex_by_gene) %in% gene_ids]
  if (length(ex_in_win)>0) {
    all_exons <- subsetByOverlaps(unlist(ex_in_win), zoom_gr)
    if (length(all_exons)>0)
      kpPlotRegions(kp, data=all_exons, col='#F5C16C', border='#8B6914',
                    r0=r0+0.01, r1=r1-0.01, lwd=0.5)
  }
  mids <- (pmax(start(win_genes), ZOOM_START) + pmin(end(win_genes), ZOOM_END))/2
  kpText(kp, chr=as.character(seqnames(win_genes)), x=mids, y=0.88,
         labels=gene_names, col='#1F3A5F', cex=0.8, font=3, r0=r0, r1=r1)
  message('  genes: ', paste(gene_names, collapse=', '))
}

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
message('Building plot...')
pdf(file.path(OUT_DIR, 'SMN2_3prime_multitrack_v5.pdf'), width=13, height=11)

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin    <- 0.18
pp$rightmargin   <- 0.04
pp$topmargin     <- 40
pp$bottommargin  <- 30
pp$ideogramheight <- 8
pp$data1inmargin  <- 4
pp$data1outmargin <- 0

kp <- plotKaryotype(genome='hg38', chromosomes=ZOOM_CHR,
                    zoom=zoom, plot.type=1, plot.params=pp,
                    main="SMN2 3' locus -- H3K9me2 / H3K27ac / WGBS methylation")
kpAddBaseNumbers(kp, tick.dist=10000, minor.tick.dist=2000,
                 add.units=TRUE, cex=0.7)

# DMR shaded band across all tracks
kpRect(kp, chr=ZOOM_CHR, x0=DMR_POS, x1=DMR_POS+300,
       y0=0, y1=1, col='#C0392B18', border=NA, r0=0, r1=1)
kpAbline(kp, v=DMR_POS, col='#C0392B', lwd=1.2, lty=2)

# Track r0/r1 layout (bottom to top):
GENE_R0 <- 0.00; GENE_R1 <- 0.07
ANN_R0  <- 0.08; ANN_R1  <- 0.17   # CpGi / H9 / cCRE (3 sub-tracks)
METH_R0 <- 0.19; METH_R1 <- 0.54   # 4 methylation lines
H27_R0  <- 0.56; H27_R1  <- 0.67   # H3K27ac peaks
H9M_R0  <- 0.69; H9M_R1  <- 1.00   # H3K9me2 bigWigs (2 conditions)

# ---- Gene track ----
draw_genes(kp, GENE_R0, GENE_R1)
kpAddLabels(kp, 'Genes', r0=GENE_R0, r1=GENE_R1, cex=0.75, col='grey30')

# ---- Annotation tracks ----
if (length(cgi)>0)
  kpPlotRegions(kp, data=cgi, col='#A8D5E2', border=NA,
                r0=ANN_R0, r1=ANN_R0+0.025)
kpAddLabels(kp, 'CpG isl', r0=ANN_R0, r1=ANN_R0+0.025, cex=0.65, col='#3A7D94')

if (length(enh)>0)
  kpPlotRegions(kp, data=enh, col='#E69F00', border=NA,
                r0=ANN_R0+0.03, r1=ANN_R0+0.055)
kpAddLabels(kp, 'H9 enh', r0=ANN_R0+0.03, r1=ANN_R0+0.055, cex=0.65, col='#E69F00')

if (length(ccre_gr)>0)
  kpPlotRegions(kp, data=ccre_gr, col='#56B4E9', border=NA,
                r0=ANN_R0+0.06, r1=ANN_R0+0.085)
kpAddLabels(kp, 'cCRE', r0=ANN_R0+0.06, r1=ANN_R0+0.085, cex=0.65, col='#56B4E9')

# ---- Methylation tracks ----
band_h <- (METH_R1 - METH_R0) / length(CONDITIONS)
for (i in seq_along(CONDITIONS)) {
  cond <- CONDITIONS[i]
  r0 <- METH_R0 + (i-1)*band_h
  r1 <- r0 + band_h - 0.005
  gr <- meth[[cond]]
  if (!is.null(gr) && length(gr)>0) {
    kpPoints(kp, data=gr, y=gr$beta,
             col=adjustcolor(PAL[cond], 0.25), pch=16, cex=0.12,
             r0=r0, r1=r1, ymin=0, ymax=1)
    sm <- rollmean_vec(gr$beta, k=15)
    ok <- !is.na(sm)
    if (sum(ok)>2) {
      gr_sm <- GRanges(ZOOM_CHR, IRanges(start(gr)[ok], width=1), beta=sm[ok])
      kpLines(kp, data=gr_sm, y=gr_sm$beta, col=PAL[cond], lwd=2,
              r0=r0, r1=r1, ymin=0, ymax=1)
    }
  }
  kpAxis(kp, ymin=0, ymax=1, r0=r0, r1=r1, side=2, cex=0.45, numticks=2)
  kpAddLabels(kp, cond, r0=r0, r1=r1, cex=0.7, col=PAL[cond])
}

# ---- H3K27ac peak tracks ----
if (length(h3k27ac_ctrl)>0)
  kpPlotRegions(kp, data=h3k27ac_ctrl, col=adjustcolor(H3K27AC_CTRL_COL,0.7),
                border=NA, r0=H27_R0, r1=H27_R0+0.045)
kpAddLabels(kp, 'H3K27ac CTRL', r0=H27_R0, r1=H27_R0+0.045,
            cex=0.7, col=H3K27AC_CTRL_COL)

if (length(h3k27ac_vpa)>0)
  kpPlotRegions(kp, data=h3k27ac_vpa, col=adjustcolor(H3K27AC_VPA_COL,0.7),
                border=NA, r0=H27_R0+0.05, r1=H27_R1)
kpAddLabels(kp, 'H3K27ac VPA', r0=H27_R0+0.05, r1=H27_R1,
            cex=0.7, col=H3K27AC_VPA_COL)

# ---- H3K9me2 bigWig tracks (smooth signal -- THE ENCODE-STYLE TRACKS) ----
h9m_tracks <- list(
  list(bw=H3K9ME2_BW$CTRL_R1, label='H3K9me2 CTR R1', col=H3K9ME2_CTRL_COL),
  list(bw=H3K9ME2_BW$CTRL_R2, label='H3K9me2 CTR R2', col=H3K9ME2_CTRL_COL),
  list(bw=H3K9ME2_BW$ASO_R1,  label='H3K9me2 ASO R1', col=H3K9ME2_ASO_COL),
  list(bw=H3K9ME2_BW$ASO_R2,  label='H3K9me2 ASO R2', col=H3K9ME2_ASO_COL)
)
for (i in seq_along(h9m_tracks)) {
  tr  <- h9m_tracks[[i]]
  if (!file.exists(tr$bw)) { message('  missing: ', tr$bw); next }
  at  <- autotrack(i, length(h9m_tracks), r0=H9M_R0, r1=H9M_R1, margin=0.08)
  kp  <- kpPlotBigWig(kp, data=tr$bw, ymax='visible.region',
                       r0=at$r0, r1=at$r1,
                       col=adjustcolor(tr$col, 0.65), border=tr$col)
  ymax <- ceiling(kp$latest.plot$computed.values$ymax * 10)/10
  kpAxis(kp, ymin=0, ymax=ymax, tick.pos=c(0, ymax),
         r0=at$r0, r1=at$r1, cex=0.45)
  kpAddLabels(kp, tr$label, r0=at$r0, r1=at$r1, cex=0.65, col=tr$col)
}
kpAddLabels(kp, 'H3K9me2', r0=H9M_R0, r1=H9M_R1,
            srt=90, pos=1, cex=0.9, label.margin=0.13, col='grey30')

# Legend
legend('topright', inset=c(0.01,0),
       legend=c('H3K9me2 CTR','H3K9me2 ASO','H3K27ac CTR','H3K27ac VPA',
                'ASO_CTRL meth','Scramble_CTRL meth','ASO_VPA meth','Scramble_VPA meth',
                'Sensitive DMR'),
       col=c(H3K9ME2_CTRL_COL, H3K9ME2_ASO_COL,
             H3K27AC_CTRL_COL, H3K27AC_VPA_COL,
             PAL['ASO_CTRL'], PAL['Scramble_CTRL'],
             PAL['ASO_VPA'], PAL['Scramble_VPA'],
             '#C0392B'),
       lty=c(1,1,1,1,1,1,1,1,2), lwd=c(2,2,2,2,2,2,2,2,1.5),
       bty='n', cex=0.6)

dev.off()
message('saved: SMN2_3prime_multitrack_v5.pdf')
