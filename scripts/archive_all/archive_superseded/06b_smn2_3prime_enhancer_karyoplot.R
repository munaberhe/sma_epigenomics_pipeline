.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(karyoploteR)
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(annotatr)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/smn2_enhancer"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# Locked SMA palette
PAL <- c(
  ASO_CTRL      = "#1F3A5F",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#C0392B",
  Scramble_VPA  = "#D4A017"
)
ENH_COL  <- "#E69F00"
CGI_COL  <- "#A8D5E2"
GENE_COL <- "#F5C16C"

# SMN2 locus - zoom to gene body + 25kb downstream (3' end focus)
SMN2_START <- 70049638
SMN2_END   <- 70078522
ZOOM_START <- SMN2_START - 5000
ZOOM_END   <- SMN2_END   + 25000
ZOOM_CHR   <- "chr5"
zoom <- toGRanges(data.frame(chr=ZOOM_CHR, start=ZOOM_START, end=ZOOM_END))
message(sprintf("Zoom: %s:%d-%d (%.1f kb)", ZOOM_CHR, ZOOM_START, ZOOM_END,
                (ZOOM_END-ZOOM_START)/1000))

# Load enhancers
message("loading enhancers...")
enh_df <- read.table(gzfile("data/reference/H9_predicted_non_promoter_non_fragments.bed.gz"),
              header=TRUE, sep="\t", stringsAsFactors=FALSE)
enh <- GRanges(enh_df$seqnames, IRanges(enh_df$start, enh_df$end))
enh <- enh[as.character(seqnames(enh))==ZOOM_CHR &
           start(enh)>=ZOOM_START & end(enh)<=ZOOM_END]
message("  ", length(enh), " enhancers in window")

# Load CpG islands (annotatr hg38)
message("loading CpG islands...")
cpg_anns <- build_annotations(genome="hg38", annotations="hg38_cpg_islands")
cgi <- cpg_anns[as.character(seqnames(cpg_anns))==ZOOM_CHR &
                start(cpg_anns)>=ZOOM_START & end(cpg_anns)<=ZOOM_END]
message("  ", length(cgi), " CpG islands in window")

# Load methylation from masked chr5 CX files
CONDITIONS <- c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA")
BY_CHR <- "results/alignments_smn1_masked/chr5_cx"

load_meth_smn2 <- function(condition) {
  files <- list.files(BY_CHR,
    pattern=sprintf("^%s_[0-9]+_chr5\.CX_report\.txt$", condition),
    full.names=TRUE)
  if (length(files)==0) return(NULL)
  message("  loading: ", condition, " (", length(files), " reps)")
  df_list <- lapply(files, function(f) {
    d <- read.table(f, header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","mC","uC","ctx","tri"))
    d <- d[d$ctx=="CG" & d$pos>=ZOOM_START & d$pos<=ZOOM_END,]
    d$beta <- d$mC/(d$mC+d$uC)
    d[d$mC+d$uC>=3, c("pos","beta")]
  })
  merged <- do.call(rbind, df_list)
  merged <- aggregate(beta~pos, merged, mean, na.rm=TRUE)
  merged <- merged[order(merged$pos),]
  GRanges(ZOOM_CHR, IRanges(merged$pos, width=1), beta=merged$beta)
}

message("loading methylation profiles...")
meth <- lapply(CONDITIONS, load_meth_smn2)
names(meth) <- CONDITIONS

# Rolling mean smoother
rollmean_vec <- function(x, k=21) {
  stats::filter(x, rep(1/k,k), sides=2)
}

# Build plot
pdf(file.path(OUT_DIR, "SMN2_3prime_enhancer_karyoplot.pdf"),
    width=12, height=9)

pp <- getDefaultPlotParams(plot.type=2)
pp$leftmargin    <- 0.20
pp$rightmargin   <- 0.04
pp$topmargin     <- 40
pp$bottommargin  <- 50
pp$ideogramheight <- 20
pp$data1height   <- 420
pp$data2height   <- 60
pp$data1inmargin <- 3
pp$data2inmargin <- 3

kp <- plotKaryotype(genome="hg38", chromosomes=ZOOM_CHR,
                    zoom=zoom, plot.type=2, plot.params=pp,
                    main="SMN2 3' enhancer locus")

kpAddBaseNumbers(kp, tick.dist=10000, tick.len=4, cex=0.6, units="Mb")

# Track positions (r0 to r1, top to bottom)
# Genes: 0.90-1.00
# CGI:   0.84-0.88
# Enh:   0.78-0.82
# Meth:  0.04-0.75 (4 conditions, each ~0.17 tall)

# Genes track - manual (makeGenesDataFromTxDb fails on small windows)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(GenomicFeatures)
all_genes <- genes(txdb)
win_genes <- all_genes[as.character(seqnames(all_genes))==ZOOM_CHR &
                        start(all_genes)<=ZOOM_END &
                        end(all_genes)>=ZOOM_START]
if (length(win_genes)>0) {
  gene_names <- tryCatch(
    mapIds(org.Hs.eg.db, keys=names(win_genes), column="SYMBOL",
           keytype="ENTREZID", multiVals="first"),
    error=function(e) setNames(names(win_genes), names(win_genes)))
  kpPlotRegions(kp, data=win_genes, col=GENE_COL, border="grey60",
                r0=0.90, r1=0.97)
  kpText(kp, data=win_genes, labels=gene_names[names(win_genes)],
         y=0.5, r0=0.90, r1=0.97, col="#1F3A5F", cex=0.7)
}
kpAddLabels(kp, "Genes", r0=0.90, r1=0.97, side="left", cex=0.75)

# CpG islands
if (length(cgi)>0)
  kpPlotRegions(kp, data=cgi, col=CGI_COL, border=NA, r0=0.84, r1=0.88)
kpAddLabels(kp, "CpG islands", r0=0.84, r1=0.88, side="left", cex=0.75)

# Enhancers
if (length(enh)>0)
  kpPlotRegions(kp, data=enh, col=ENH_COL, border=ENH_COL, r0=0.82, r1=0.86)
kpAddLabels(kp, "H9 enhancers", r0=0.82, r1=0.86, side="left", cex=0.75)

# ENCODE cCRE enhancers track (pELS + dELS)
message("loading ENCODE cCREs...")
cCRE_df <- read.table("data/reference/encode_cCREs_hg38.bed",
                       header=FALSE, sep="\t", stringsAsFactors=FALSE,
                       col.names=c("chr","start","end","id1","id2","type"))
cCRE_enh <- cCRE_df[grepl("ELS", cCRE_df$type) &
                     cCRE_df$chr==ZOOM_CHR &
                     cCRE_df$start>=ZOOM_START &
                     cCRE_df$end<=ZOOM_END,]
message("  ", nrow(cCRE_enh), " ENCODE enhancer-like elements in window")
if (nrow(cCRE_enh)>0) {
  cCRE_gr <- GRanges(cCRE_enh$chr, IRanges(cCRE_enh$start, cCRE_enh$end))
  kpPlotRegions(kp, data=cCRE_gr, col="#56B4E9", border="#56B4E9",
                r0=0.77, r1=0.81)
}
kpAddLabels(kp, "ENCODE cCREs", r0=0.77, r1=0.81, side="left", cex=0.75)

# DMR position marker
kpAbline(kp, v=70088223, col="#C0392B", lwd=1.5, lty=2)

# Methylation tracks (4 conditions)
band_h <- 0.165
for (i in seq_along(CONDITIONS)) {
  cond <- CONDITIONS[i]
  r0 <- 0.04 + (i-1)*band_h
  r1 <- r0 + band_h - 0.01
  gr <- meth[[cond]]
  if (!is.null(gr) && length(gr)>0) {
    kpPoints(kp, data=gr, y=gr$beta,
             col=adjustcolor(PAL[cond],0.3), pch=16, cex=0.15,
             r0=r0, r1=r1, ymin=0, ymax=1)
    sm <- rollmean_vec(gr$beta, k=15)
    ok <- !is.na(sm)
    if (sum(ok)>2) {
      gr_sm <- GRanges(ZOOM_CHR, IRanges(start(gr)[ok], width=1), beta=sm[ok])
      kpLines(kp, data=gr_sm, y=gr_sm$beta,
              col=PAL[cond], lwd=2,
              r0=r0, r1=r1, ymin=0, ymax=1)
    }
  }
  kpAxis(kp, ymin=0, ymax=1, r0=r0, r1=r1, side=2,
         cex=0.5, numticks=3)
  kpAddLabels(kp, cond, r0=r0, r1=r1, side="left", cex=0.65)
}

# DMR track (data panel 2, below ideogram)
dmr_csv <- "results/smn2_enhancer/SMN2_region_DMR_enhancer_overlaps.csv"
if (file.exists(dmr_csv)) {
  dmrs <- read.csv(dmr_csv)
  dmrs <- dmrs[dmrs$seqnames==ZOOM_CHR &
               dmrs$start>=ZOOM_START & dmrs$end<=ZOOM_END,]
  if (nrow(dmrs)>0) {
    dmr_gr <- GRanges(dmrs$seqnames, IRanges(dmrs$start, dmrs$end))
    dmr_col <- ifelse(dmrs$methDiff>0, "#C0392B", "#1F3A5F")
    kpPlotRegions(kp, data=dmr_gr, col=dmr_col, border=NA,
                  data.panel=2, r0=0.1, r1=0.9)
    kpAddLabels(kp, "DMRs", data.panel=2, r0=0.1, r1=0.9,
                side="left", cex=0.7)
  }
}

dev.off()
message("saved: SMN2_3prime_enhancer_karyoplot.pdf")
