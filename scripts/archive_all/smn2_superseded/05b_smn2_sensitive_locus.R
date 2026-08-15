.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
source("scripts/pipeline/00_sma_palette.R")
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

CHR5_MASKED <- "results/alignments_smn1_masked/chr5_cx"
SENS_CSV    <- "results/smn2_local_dmr/SMN2_sensitive_DMRs_all_contrasts.csv"
OUT_DIR     <- "results/smn2_sensitive_locus"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FLANK    <- 5000
WIN_SIZE <- 300

LOCI <- list(SMN2=list(chr="chr5", start=70049638, end=70078522, strand="+"))
EXONS <- list(SMN2=data.frame(
  exon=c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
  start=c(70049638,70063415,70065965,70066934,70067294,70069235,70070641,70076521,70077019),
  end  =c(70049766,70063486,70066084,70067134,70067446,70069330,70070751,70076574,70077592),
  is_target=c(F,F,F,F,F,F,F,T,F)))

COMPARISONS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",    cond1="ASO_CTRL",    cond2="Scramble_CTRL"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",cond1="Scramble_VPA",cond2="Scramble_CTRL"),
  list(name="ASO_VPA_vs_Scramble_CTRL",     cond1="ASO_VPA",     cond2="Scramble_CTRL"),
  list(name="ASO_VPA_vs_ASO_CTRL",          cond1="ASO_VPA",     cond2="ASO_CTRL"),
  list(name="ASO_VPA_vs_Scramble_VPA",      cond1="ASO_VPA",     cond2="Scramble_VPA")
)
NEEDED <- unique(unlist(lapply(COMPARISONS, function(x) c(x$cond1, x$cond2))))
COND_COLOURS <- c(ASO_CTRL="#1F3A5F", Scramble_CTRL="#6B7280",
                  ASO_VPA="#C0392B",  Scramble_VPA="#D4A017")

# Load sensitive DMRs
if (!file.exists(SENS_CSV)) stop("Run 07_smn2_sensitive.R first: ", SENS_CSV)
sdf <- read.csv(SENS_CSV, stringsAsFactors=FALSE)
SENS_DMRS <- list()
for (ct in unique(sdf$contrast)) {
  sub <- sdf[sdf$contrast==ct,]
  # keep only the most significant hit per genomic locus (by start)
  sub <- sub[order(sub$pValue),]
  sub <- sub[!duplicated(sub$start),]
  if (nrow(sub)>0)
    SENS_DMRS[[ct]] <- GRanges(seqnames=sub$seqnames,
      ranges=IRanges(sub$start, sub$end),
      regionType=sub$regionType, pValue=sub$pValue)
}
message("Loaded sensitive DMRs for ", length(SENS_DMRS), " contrasts")

build_gff <- function() {
  rows <- list()
  l <- LOCI[["SMN2"]]
  rows[[1]] <- data.frame(chr=l$chr, start=l$start, end=l$end,
    strand=l$strand, type="gene", name="SMN2", stringsAsFactors=FALSE)
  ex <- EXONS[["SMN2"]]
  for (i in seq_len(nrow(ex)))
    rows[[length(rows)+1]] <- data.frame(chr=l$chr, start=ex$start[i], end=ex$end[i],
      strand=l$strand, type="exon", name=sprintf("SMN2_%s",ex$exon[i]), stringsAsFactors=FALSE)
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start,df$end),
          strand=df$strand, type=df$type, name=df$name)
}
GEs <- build_gff()

read_masked_cpg <- function(condition) {
  message("  loading: ", condition)
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

message("Loading pooled CpG data...")
pooled <- lapply(NEEDED, read_masked_cpg); names(pooled) <- NEEDED

for (ct in COMPARISONS) {
  dmrs <- SENS_DMRS[[ct$name]]
  n_dmr <- if(is.null(dmrs)) 0 else length(dmrs)
  fname <- sprintf("SMN2_sensitive_%s.pdf", ct$name)
  message("  plotting: ", fname, " (", n_dmr, " DMRs)")
  cairo_pdf(file.path(OUT_DIR, fname), width=12, height=7)
  par(mar=c(5,4,3,1)+0.5, cex=0.9, bg="white",
      col.axis="black", col.lab="black", col.main="black", fg="black")
  locus  <- LOCI[["SMN2"]]
  region <- GRanges(seqnames=locus$chr,
                    ranges=IRanges(locus$start-FLANK, locus$end+FLANK))
  plotLocalMethylationProfile(
    methylationData1=pooled[[ct$cond1]],
    methylationData2=pooled[[ct$cond2]],
    region=region,
    DMRs=if(n_dmr>0) list("sensitive DMRs"=dmrs) else NULL,
    conditionsNames=c(ct$cond1, ct$cond2),
    gff=GEs, windowSize=WIN_SIZE, context="CG",
    col=NULL,
    main=sprintf("SMN2 sensitive scan: %s vs %s (%d DMRs)",
                 ct$cond1, ct$cond2, n_dmr),
    plotMeanLines=TRUE, plotPoints=TRUE)
  ex <- EXONS[["SMN2"]]
  for (i in seq_len(nrow(ex)))
    mtext(ex$exon[i], side=1, at=(ex$start[i]+ex$end[i])/2, line=0.8, cex=0.5,
          col=if(ex$is_target[i])"red" else "black",
          font=if(ex$is_target[i])2 else 1)
  dev.off()
}
message("\nDone. Outputs in: ", OUT_DIR)
