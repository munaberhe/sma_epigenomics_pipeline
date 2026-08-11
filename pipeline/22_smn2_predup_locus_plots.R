#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

CX_DIR   <- "results/smn2_predup_cx"
OUT_SMN2 <- "results/smn2_locus_final"
OUT_SENS <- "results/smn2_sensitive_locus"
dir.create(OUT_SMN2, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_SENS, recursive=TRUE, showWarnings=FALSE)

SMN2_REG <- GRanges("chr5", IRanges(70040000, 70090000))
WIN_SIZE <- 300
DMR_COL  <- "#E31A1C"

COND_COLS <- c(
  ASO_CTRL      = "#1F3A5F",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#C0392B",
  Scramble_VPA  = "#F0A500"
)

CONDITIONS <- c("ASO_CTRL","ASO_VPA","Scramble_CTRL","Scramble_VPA")

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond_a="Scramble_CTRL", cond_b="ASO_CTRL",
       label="ASO effect (CTRL background)"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       cond_a="Scramble_CTRL", cond_b="Scramble_VPA",
       label="VPA effect (HDAC inhibitor)"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       cond_a="Scramble_VPA", cond_b="ASO_VPA",
       label="ASO effect on VPA background"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       cond_a="ASO_CTRL", cond_b="ASO_VPA",
       label="VPA effect on ASO background")
)

LOCI_COORDS <- list(
  SMN1 = list(chr="chr5", start=70924941, end=70953015, strand="+"),
  SMN2 = list(chr="chr5", start=70049638, end=70078522, strand="+")
)
EXONS <- list(
  SMN1 = data.frame(
    exon=c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start=c(70925030,70938807,70941357,70942326,70942686,
            70944627,70946033,70951913,70952411),
    end=c(70925158,70938878,70941476,70942526,70942838,
          70944722,70946143,70951966,70952984),
    is_target=c(F,F,F,F,F,F,F,T,F)
  ),
  SMN2 = data.frame(
    exon=c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start=c(70049638,70063415,70065965,70066934,70067294,
            70069235,70070641,70076521,70077019),
    end=c(70049766,70063486,70066084,70067134,70067446,
          70069330,70070751,70076574,70077592),
    is_target=c(F,F,F,F,F,F,F,T,F)
  )
)


message("Loading pre-dedup SMN2 data...")
load_smn2 <- function(condition) {
  message("  ", condition)
  grs <- lapply(1:3, function(r) {
    pattern <- paste0(condition, "_", r, "_*_smn2.txt")
    f <- list.files(CX_DIR, pattern=glob2rx(pattern), full.names=TRUE)
    if (length(f) == 0) { message("    missing"); return(NULL) }
    d <- read.table(f[1], header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","countM","countU","context","tri"),
      colClasses=c("character","integer","character","integer",
                   "integer","character","character"))
    if (nrow(d) == 0) return(NULL)
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos, d$pos),
            strand=d$strand, readsM=d$countM,
            readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  grs <- Filter(Negate(is.null), grs)
  if (length(grs) == 0) return(NULL)
  poolMethylationDatasets(GRangesList(grs))
}

meth_predup <- lapply(CONDITIONS, load_smn2)
names(meth_predup) <- CONDITIONS

for (cond in CONDITIONS) {
  m <- meth_predup[[cond]]
  if (is.null(m)) next
  message(cond, ": ", length(m), " CpGs | >=3x: ", sum(m$readsN>=3),
          " | >=5x: ", sum(m$readsN>=5), " | median: ", round(median(m$readsN),1))
}


SENSITIVE_DMRS <- list()
sens_csv <- "results/smn2_local_dmr/SMN2_sensitive_DMRs_all_contrasts.csv"
if (file.exists(sens_csv)) {
  sdf <- read.csv(sens_csv)
  for (ct in unique(sdf$contrast)) {
    sub <- sdf[sdf$contrast==ct,]
    if (nrow(sub)>0)
      SENSITIVE_DMRS[[ct]] <- GRanges(sub$seqnames,
        IRanges(sub$start,sub$end),
        regionType=sub$regionType, pValue=sub$pValue)
  }
}

# Also load pre-dedup canonical DMRs
predup_dmrs <- readRDS("results/smn2_predup_comparison/dmr_canonical_predup_SMN2.rds")


build_smn_gff <- function() {
  rows <- list()
  for (ln in names(LOCI_COORDS)) {
    l <- LOCI_COORDS[[ln]]
    rows[[length(rows)+1]] <- data.frame(
      chr=l$chr, start=l$start, end=l$end,
      strand=l$strand, type="gene", name=ln)
    ex <- EXONS[[ln]]
    for (i in seq_len(nrow(ex)))
      rows[[length(rows)+1]] <- data.frame(
        chr=l$chr, start=ex$start[i], end=ex$end[i],
        strand=l$strand, type="exon",
        name=sprintf("%s_%s", ln, ex$exon[i]))
  }
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start,df$end),
          strand=df$strand, type=df$type, name=df$name)
}
SMN_GFF <- build_smn_gff()


plot_smn2 <- function(meth_a, meth_b, ct, dmrs=NULL,
                      annotation=NULL, title=NULL) {
  region <- GRanges("chr5", IRanges(70049638-2000, 70078522+2000))
  ma <- meth_a[!is.na(meth_a$readsN) & meth_a$readsN >= 3]
  mb <- meth_b[!is.na(meth_b$readsN) & meth_b$readsN >= 3]

  plotLocalMethylationProfile(
    methylationData1 = ma,
    methylationData2 = mb,
    region           = region,
    DMRs             = NULL,
    conditionsNames  = c(ct$cond_a, ct$cond_b),
    gff              = SMN_GFF,
    windowSize       = WIN_SIZE,
    context          = "CG",
    col              = c(COND_COLS[ct$cond_a], COND_COLS[ct$cond_b]),
    main             = if (!is.null(title)) title else
                       sprintf("SMN2 (pre-dedup): %s vs %s\n%s",
                               ct$cond_a, ct$cond_b, ct$label),
    plotMeanLines    = TRUE,
    plotPoints       = TRUE  # show all CpG dots
  )

  # DMR overdraw
  if (!is.null(dmrs) && length(dmrs) > 0 && length(dmrs) <= 30) {
    tryCatch(par(mfg=c(1,1)), error=function(e) NULL)
    usr <- par("usr")
    meth_top <- usr[4] - (usr[4]-usr[3]) * 0.03
    meth_bot <- usr[4] - (usr[4]-usr[3]) * 0.10
    dmr_df <- as.data.frame(dmrs)
    for (i in seq_len(nrow(dmr_df))) {
      rect(xleft=dmr_df$start[i]-500, xright=dmr_df$end[i]+500,
           ybottom=meth_bot, ytop=meth_top,
           col=DMR_COL, border=DMR_COL, lwd=0.5)
    }
  }

  # Exon labels
  ex <- EXONS[["SMN2"]]
  for (i in seq_len(nrow(ex)))
    mtext(ex$exon[i], side=1, at=(ex$start[i]+ex$end[i])/2,
          line=-1.5, cex=0.45,
          col=if(ex$is_target[i]) DMR_COL else "black",
          font=if(ex$is_target[i]) 2 else 1)

  if (!is.null(annotation))
    mtext(annotation, side=3, line=0.2, cex=0.7, col="grey40", adj=1)
}


message("\n=== SMN2 locus plots (pre-dedup) ===")
for (ct in CONTRASTS) {
  fname <- paste0("SMN_locus_predup_", ct$name, ".pdf")
  message("  ", fname)
  cairo_pdf(file.path(OUT_SMN2, fname), width=11, height=7,
            bg="white", onefile=TRUE)
  par(mfrow=c(1,1), bg="white", col.axis="black",
      col.lab="black", col.main="black", fg="black")
  dmrs <- predup_dmrs[[ct$name]]
  tryCatch(
    plot_smn2(meth_predup[[ct$cond_a]], meth_predup[[ct$cond_b]],
              ct, dmrs=dmrs,
              annotation="Pre-dedup | plotPoints=TRUE | min coverage=3x"),
    error=function(e) { plot.new(); text(0.5,0.5,e$message,col="red") }
  )
  dev.off()
}


message("\n=== SMN2 sensitive plots (pre-dedup) ===")
for (ct in CONTRASTS) {
  fname <- paste0("SMN2_sensitive_predup_", ct$name, ".pdf")
  message("  ", fname)
  cairo_pdf(file.path(OUT_SENS, fname), width=11, height=7,
            bg="white", onefile=TRUE)
  par(mfrow=c(1,1), bg="white", col.axis="black",
      col.lab="black", col.main="black", fg="black")
  sens_dmrs <- SENSITIVE_DMRS[[ct$name]]
  tryCatch(
    plot_smn2(meth_predup[[ct$cond_a]], meth_predup[[ct$cond_b]],
              ct, dmrs=sens_dmrs,
              title=sprintf("SMN2 sensitive (pre-dedup): %s vs %s\n%s",
                            ct$cond_a, ct$cond_b, ct$label),
              annotation="Pre-dedup | sensitive params | min coverage=3x"),
    error=function(e) { plot.new(); text(0.5,0.5,e$message,col="red") }
  )
  dev.off()
}

message("\n=== ALL DONE ===")
