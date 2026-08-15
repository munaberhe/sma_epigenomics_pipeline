suppressPackageStartupMessages({
  library(karyoploteR); library(GenomicRanges)
  library(rtracklayer); library(DMRcaller)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
BY_CHR_DIR <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/smn2_enhancer"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

COND_COLOURS <- c(Scramble_CTRL="#6B7280", ASO_CTRL="#1F3A5F",
                  Scramble_VPA="#D4A017", ASO_VPA="#C0392B")
COND_ORDER <- c("Scramble_CTRL","ASO_CTRL","Scramble_VPA","ASO_VPA")
REPS <- 1:3
LOCUS_CHR <- "chr5"; LOCUS_START <- 70083000; LOCUS_END <- 70095000
SMN2_START <- 70049638; SMN2_END <- 70077595
DMR_START <- 70088223; DMR_END <- 70088522; DMR_DELTA <- 36.3

# Load enhancers
load_bed_window <- function(path, header, col_names=NULL) {
  if (!file.exists(path)) return(GRanges())
  d <- if (header) read.table(gzfile(path), header=TRUE, sep="\t", stringsAsFactors=FALSE) else {
    x <- read.table(gzfile(path), header=FALSE, sep="\t", stringsAsFactors=FALSE)
    if (!is.null(col_names)) colnames(x) <- col_names; x }
  chrcol <- if ("seqnames" %in% colnames(d)) "seqnames" else "chr"
  keep <- d[[chrcol]]==LOCUS_CHR & d$end>=LOCUS_START & d$start<=LOCUS_END
  d <- d[keep,,drop=FALSE]
  if (nrow(d)==0) return(GRanges())
  GRanges(d[[chrcol]], IRanges(d$start, d$end))
}

h9_gr <- load_bed_window("data/reference/H9_predicted_non_promoter_non_fragments.bed.gz", TRUE)

# E081 fetal brain - not available (network restricted); set empty
e081_gr <- GRanges()
message("E081 fetal brain track: not available on HPC network")

# ENCODE cCREs
ccre_df <- read.table("data/reference/encode_cCREs_hg38.bed", header=FALSE, sep="\t",
                       col.names=c("chr","start","end","id1","id2","type"))
ccre_sub <- ccre_df[grepl("ELS",ccre_df$type) & ccre_df$chr==LOCUS_CHR &
                    ccre_df$start>=LOCUS_START & ccre_df$end<=LOCUS_END,]
ccre_gr <- if (nrow(ccre_sub)>0) GRanges(ccre_sub$chr, IRanges(ccre_sub$start, ccre_sub$end)) else GRanges()

message(sprintf("H9 enhancers: %d | E081 fetal brain: %d | ENCODE cCREs: %d",
                length(h9_gr), length(e081_gr), length(ccre_gr)))

# Methylation profiles
read_cpg <- function(f) {
  d <- read.table(gzfile(f), header=FALSE, sep="\t",
                  col.names=c("chr","pos","strand","countM","countU","context","tri"))
  d <- d[d$context=="CG",]
  GRanges(d$chr, IRanges(d$pos,d$pos), strand=d$strand,
          readsM=d$countM, readsN=d$countM+d$countU,
          context=d$context, trinucleotide_context=d$tri)
}
message("Loading methylation profiles...")
profiles <- lapply(COND_ORDER, function(cond) {
  files <- file.path(BY_CHR_DIR, sprintf("%s_%d_chr5.CpG_report.txt.gz", cond, REPS))
  files <- files[file.exists(files)]
  if (length(files)==0) return(NULL)
  pooled <- poolMethylationDatasets(GRangesList(lapply(files, read_cpg)))
  region <- GRanges(LOCUS_CHR, IRanges(LOCUS_START, LOCUS_END))
  prof <- computeMethylationProfile(pooled, region, windowSize=300, context="CG")
  prof[prof$readsN>=2]
})
names(profiles) <- COND_ORDER

# Plot
pdf(file.path(OUT_DIR, "SMN2_3prime_enhancer_karyoplot_v3.pdf"), width=13, height=9)
zoom_gr <- toGRanges(data.frame(chr=LOCUS_CHR, start=LOCUS_START, end=LOCUS_END))
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin=0.16; pp$rightmargin=0.06; pp$topmargin=30
pp$bottommargin=30; pp$ideogramheight=8; pp$data1inmargin=6

kp <- plotKaryotype(genome="hg38", chromosomes=LOCUS_CHR, zoom=zoom_gr,
                    plot.params=pp, cex=0.9)
kpAddBaseNumbers(kp, tick.dist=10000, tick.len=5, cex=0.7,
                 minor.tick.dist=2000, minor.tick.len=2, add.units=TRUE)
kpAddMainTitle(kp, main="SMN2 3' downstream — sensitive-reanalysis DMR (+36.3%, ASO_VPA vs Scramble_CTRL)", cex=1.0)

# DMR shaded band
kpRect(kp, chr=LOCUS_CHR, x0=DMR_START, x1=DMR_END, y0=0, y1=1,
       col="#C0392B22", border=NA, r0=0, r1=1)

# Annotation tracks
ann_r0 <- 0.78; ann_r1 <- 0.98
n_ann <- 5

at <- autotrack(1, n_ann, r0=ann_r0, r1=ann_r1, margin=0.15)
kpRect(kp, chr=LOCUS_CHR, x0=SMN2_START, x1=SMN2_END, y0=0.1, y1=0.9,
       col="#E8B142", border="#9C7100", r0=at$r0, r1=at$r1, lwd=0.8)
kpText(kp, chr=LOCUS_CHR, x=(SMN2_START+SMN2_END)/2, y=0.5, labels="SMN2",
       cex=0.85, font=3, r0=at$r0, r1=at$r1, col="#3D2C00")
kpAddLabels(kp, "Genes", r0=at$r0, r1=at$r1, cex=0.8)

at <- autotrack(2, n_ann, r0=ann_r0, r1=ann_r1, margin=0.15)
kpAddLabels(kp, "H9 enh (ESC)", r0=at$r0, r1=at$r1, cex=0.75)
if (length(h9_gr)>0) kpRect(kp, data=h9_gr, y0=0.25, y1=0.75, col="#5FB6E5", border=NA, r0=at$r0, r1=at$r1) else
  kpText(kp, chr=LOCUS_CHR, x=(LOCUS_START+LOCUS_END)/2, y=0.5, labels="(none)", cex=0.6, col="#888", r0=at$r0, r1=at$r1)

at <- autotrack(3, n_ann, r0=ann_r0, r1=ann_r1, margin=0.15)
kpAddLabels(kp, "Fetal brain (E081)", r0=at$r0, r1=at$r1, cex=0.75)
if (length(e081_gr)>0) kpRect(kp, data=e081_gr, y0=0.25, y1=0.75, col="#4E944F", border=NA, r0=at$r0, r1=at$r1) else
  kpText(kp, chr=LOCUS_CHR, x=(LOCUS_START+LOCUS_END)/2, y=0.5, labels="(none)", cex=0.6, col="#888", r0=at$r0, r1=at$r1)

at <- autotrack(4, n_ann, r0=ann_r0, r1=ann_r1, margin=0.15)
kpAddLabels(kp, "ENCODE cCREs", r0=at$r0, r1=at$r1, cex=0.75)
if (length(ccre_gr)>0) kpRect(kp, data=ccre_gr, y0=0.25, y1=0.75, col="#A6589A", border=NA, r0=at$r0, r1=at$r1) else
  kpText(kp, chr=LOCUS_CHR, x=(LOCUS_START+LOCUS_END)/2, y=0.5, labels="(none)", cex=0.6, col="#888", r0=at$r0, r1=at$r1)

# Methylation tracks
meth_r0 <- 0.05; meth_r1 <- 0.72
for (i in seq_along(COND_ORDER)) {
  cond <- COND_ORDER[i]
  at <- autotrack(i, length(COND_ORDER), r0=meth_r0, r1=meth_r1, margin=0.18)
  prof <- profiles[[cond]]
  if (is.null(prof) || length(prof)==0) next
  mid <- (start(prof)+end(prof))/2
  kpAddLabels(kp, cond, r0=at$r0, r1=at$r1, cex=0.85)
  kpAxis(kp, r0=at$r0, r1=at$r1, ymin=0, ymax=1, numticks=3, cex=0.55)
  kpLines(kp, chr=LOCUS_CHR, x=mid, y=prof$Proportion,
          r0=at$r0, r1=at$r1, ymin=0, ymax=1,
          col=COND_COLOURS[cond], lwd=1.7)
}

kpText(kp, chr=LOCUS_CHR, x=(DMR_START+DMR_END)/2, y=0.95,
       labels=sprintf("DMR +%.1f%%", DMR_DELTA), cex=0.8, col="#7B1F12", font=2,
       r0=ann_r0, r1=ann_r1)
dev.off()
message("saved: SMN2_3prime_enhancer_karyoplot_v3.pdf")
