.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})
source("scripts/pipeline/00_sma_palette.R")
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

CHR5_MASKED   <- "results/alignments_smn1_masked/chr5_cx"
BY_CHR_UNMASK <- "results/alignments/bs/by_chr"
OUT_DIR       <- "results/smn2_locus_final"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FLANK    <- 2000
WIN_SIZE <- 300

LOCI <- list(
  SMN1 = list(chr="chr5", start=70924941, end=70953015, strand="+"),
  SMN2 = list(chr="chr5", start=70049638, end=70078522, strand="+")
)
EXONS <- list(
  SMN1 = data.frame(
    exon=c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start=c(70925030,70938807,70941357,70942326,70942686,70944627,70946033,70951913,70952411),
    end  =c(70925158,70938878,70941476,70942526,70942838,70944722,70946143,70951966,70952984),
    is_target=c(F,F,F,F,F,F,F,T,F)),
  SMN2 = data.frame(
    exon=c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start=c(70049638,70063415,70065965,70066934,70067294,70069235,70070641,70076521,70077019),
    end  =c(70049766,70063486,70066084,70067134,70067446,70069330,70070751,70076574,70077592),
    is_target=c(F,F,F,F,F,F,F,T,F))
)
COMPARISONS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",       cond1="ASO_CTRL",    cond2="Scramble_CTRL"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",    cond1="Scramble_VPA",cond2="Scramble_CTRL"),
  list(name="ASO_VPA_vs_Scramble_CTRL",         cond1="ASO_VPA",     cond2="Scramble_CTRL"),
  list(name="ASO_VPA_vs_ASO_CTRL",              cond1="ASO_VPA",     cond2="ASO_CTRL"),
  list(name="ASO_VPA_vs_Scramble_VPA",          cond1="ASO_VPA",     cond2="Scramble_VPA")
)
NEEDED <- unique(unlist(lapply(COMPARISONS, function(x) c(x$cond1, x$cond2))))
COND_COLOURS <- c(ASO_CTRL="#1F3A5F", Scramble_CTRL="#6B7280",
                  ASO_VPA="#C0392B",  Scramble_VPA="#D4A017")

build_gff <- function() {
  rows <- list()
  for (ln in names(LOCI)) {
    l <- LOCI[[ln]]
    rows[[length(rows)+1]] <- data.frame(chr=l$chr, start=l$start, end=l$end,
      strand=l$strand, type="gene", name=ln, stringsAsFactors=FALSE)
    ex <- EXONS[[ln]]
    for (i in seq_len(nrow(ex)))
      rows[[length(rows)+1]] <- data.frame(chr=l$chr, start=ex$start[i], end=ex$end[i],
        strand=l$strand, type="exon", name=sprintf("%s_%s",ln,ex$exon[i]), stringsAsFactors=FALSE)
  }
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start,df$end),
          strand=df$strand, type=df$type, name=df$name)
}
GEs <- build_gff()

read_masked_cpg <- function(condition) {
  message("  masked: ", condition)
  grs <- lapply(1:3, function(r) {
    path <- file.path(CHR5_MASKED, paste0(condition,"_",r,"_chr5.CX_report.txt"))
    d <- read.table(path, header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","countM","countU","context","tri"),
      colClasses=c("character","integer","character","integer","integer","character","character"))
    d <- d[d$context=="CG",]
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos,d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

read_unmasked_cpg <- function(condition) {
  message("  unmasked: ", condition)
  files <- file.path(BY_CHR_UNMASK, sprintf("%s_%d_chr5.CpG_report.txt.gz",condition,1:3))
  files <- files[file.exists(files)]
  grs <- lapply(files, function(f) {
    d <- read.table(gzfile(f), header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","countM","countU","context","tri"),
      colClasses=c("character","integer","character","integer","integer","character","character"))
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos,d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

# CLEAN plot - NO DMR overlay, NO sensitive DMR title
plot_one <- function(pooled, comp, locus_name) {
  locus  <- LOCI[[locus_name]]
  region <- GRanges(seqnames=locus$chr,
                    ranges=IRanges(locus$start-FLANK, locus$end+FLANK))
  m1 <- pooled[[comp$cond1]]; 
  m2 <- pooled[[comp$cond2]]; 
  plotLocalMethylationProfile(
    methylationData1=m1, methylationData2=m2,
    region=region, DMRs=NULL,
    conditionsNames=c(comp$cond1, comp$cond2),
    gff=GEs, windowSize=WIN_SIZE, context="CG",
    col=NULL,
    main=sprintf("%s: %s vs %s", locus_name, comp$cond1, comp$cond2),
    plotMeanLines=TRUE, plotPoints=TRUE)
  ex <- EXONS[[locus_name]]
  for (i in seq_len(nrow(ex)))
    mtext(ex$exon[i], side=1, at=(ex$start[i]+ex$end[i])/2, line=-1.5, cex=0.45,
          col=if(ex$is_target[i])"red" else "black",
          font=if(ex$is_target[i])2 else 1)
  usr <- par("usr")
  text(usr[1], usr[3]+(usr[4]-usr[3])*0.05, labels=locus_name, cex=0.8, font=2, adj=c(0,0.5))
}

for (alignment in c("masked","unmasked")) {
  message("\n", alignment, " alignment")
  reader <- if(alignment=="masked") read_masked_cpg else read_unmasked_cpg
  pooled <- lapply(NEEDED, reader); names(pooled) <- NEEDED
  for (ct in COMPARISONS) {
    fname <- sprintf("SMN_locus_%s_%s.pdf", alignment, ct$name)
    message("  plotting: ", fname)
    h  <- if(alignment=="masked") 7 else 9.5
    nr <- if(alignment=="masked") 1 else 2
    cairo_pdf(file.path(OUT_DIR, fname), width=12, height=h)
    par(mfrow=c(nr,1), mar=c(5,4,3,1)+0.5, cex=0.9,
        bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")
    if(alignment=="unmasked") plot_one(pooled, ct, "SMN1")
    plot_one(pooled, ct, "SMN2")
    dev.off()
  }
  fname_all <- sprintf("SMN_locus_%s_all_comparisons.pdf", alignment)
  message("  plotting combined: ", fname_all)
  n_loci   <- if(alignment=="masked") 1 else 2
  n_panels <- length(COMPARISONS)*n_loci
  n_rows   <- ceiling(n_panels/2)
  cairo_pdf(file.path(OUT_DIR, fname_all), width=16, height=n_rows*4.5)
  par(mfrow=c(n_rows,2), mar=c(5,4,3,1)+0.5, cex=0.75,
      bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")
  for (ct in COMPARISONS) {
    if(alignment=="unmasked") plot_one(pooled, ct, "SMN1")
    plot_one(pooled, ct, "SMN2")
  }
  dev.off()
}
message("\nDone. Outputs in: ", OUT_DIR)
