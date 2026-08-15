.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(data.table)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

# ---- Answering Alberto's question directly ----
# "Did you check if there's an enhancer in one of the introns that is now
# contacting the promoter and increasing activation of the gene? Whether in
# the last 2 introns, are there any enhancers?"
#
# SMN2 transcript ENST00000380743.9 (full-length, includes exon 7):
#   Exon 7 (ASO target): chr5:70,070,641-70,070,751
#   Intron 7 (last-but-one): chr5:70,070,752-70,076,520 (~5.8kb)
#   Intron 8 (last): chr5:70,076,575-70,077,018 (~450bp)
#   Combined window with small flank: chr5:70,069,752-70,078,018
#
# Checks: H9 ESC predicted enhancers, ENCODE candidate cis-regulatory
# elements (cCREs), CpG islands -- against this exact window, plus
# four-condition methylation across the same region.

BY_CHR_UNMASK <- "results/alignments/bs/by_chr"
OUT_DIR       <- "results/smn2_enhancer"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

EXON7_START <- 70070641; EXON7_END <- 70070751
INTRON7_START <- 70070752; INTRON7_END <- 70076520
INTRON8_START <- 70076575; INTRON8_END <- 70077018
EXON8_START <- 70076521; EXON8_END <- 70076574
EXON9_START <- 70077019; EXON9_END <- 70077595

FLANK <- 1000
PLOT_START <- EXON7_START - FLANK
PLOT_END   <- EXON9_END + FLANK
ZOOM_CHR   <- "chr5"

COND_COLOURS <- c(
  ASO_CTRL      = "#1B4F8A",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500"
)

# ---- Annotation checks: H9 enhancers, ENCODE cCREs, CpG islands ----
message("Checking H9 predicted enhancers in introns 7-8 window...")
enh_df <- fread("data/reference/H9_predicted_non_promoter_non_fragments.bed.gz", header=TRUE)
enh_gr <- GRanges(enh_df$seqnames, IRanges(enh_df$start, enh_df$end))
window_gr <- GRanges(ZOOM_CHR, IRanges(PLOT_START, PLOT_END))
enh_in_window <- subsetByOverlaps(enh_gr, window_gr)
message("  H9 enhancers overlapping intron 7-8 window: ", length(enh_in_window))
if (length(enh_in_window) > 0) print(enh_in_window)

message("\nChecking ENCODE cCREs in introns 7-8 window...")
ccre_df <- read.table("data/reference/encode_cCREs_hg38.bed", header=FALSE, sep="\t",
                       col.names=c("chr","start","end","id1","id2","type"))
ccre_sub <- ccre_df[ccre_df$chr==ZOOM_CHR & ccre_df$end>=PLOT_START & ccre_df$start<=PLOT_END, ]
message("  All ENCODE cCREs (any type) in window: ", nrow(ccre_sub))
if (nrow(ccre_sub) > 0) print(ccre_sub)
ccre_enh <- ccre_sub[grepl("ELS", ccre_sub$type), ]
message("  ENCODE enhancer-like signatures (ELS) specifically: ", nrow(ccre_enh))
ccre_gr <- if (nrow(ccre_sub)>0) GRanges(ccre_sub$chr, IRanges(ccre_sub$start, ccre_sub$end)) else GRanges()

message("\n=== SUMMARY: does the public annotation support an enhancer in SMN2 introns 7-8? ===")
if (length(enh_in_window)==0 && nrow(ccre_enh)==0) {
  message("No H9 ESC predicted enhancer and no ENCODE ELS signature found in this window.")
  message("This does not rule out a primed/cell-type-specific enhancer not captured by")
  message("standard annotation models -- consistent with Alberto's point that his lab's")
  message("models find more enhancers than standard methods.")
} else {
  message("At least one candidate enhancer signature found -- see details above.")
}

# ---- Methylation across the window, all four conditions ----
CONDITIONS <- c("ASO_CTRL", "Scramble_CTRL", "ASO_VPA", "Scramble_VPA")
read_unmasked_cpg <- function(condition) {
  message("  loading: ", condition)
  files <- file.path(BY_CHR_UNMASK,
                     sprintf("%s_%d_chr5.CpG_report.txt.gz", condition, 1:3))
  files <- files[file.exists(files)]
  if (length(files) == 0) return(NULL)
  grs <- lapply(files, function(f) {
    d <- read.table(gzfile(f), header=FALSE, sep="\t",
                    col.names=c("chr","pos","strand","countM","countU","context","tri"),
                    colClasses=c("character","integer","character","integer",
                                 "integer","character","character"))
    d <- d[d$context=="CG", ]
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos,d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

message("\nLoading chr5 methylation for all 4 conditions...")
pooled <- setNames(lapply(CONDITIONS, read_unmasked_cpg), CONDITIONS)

# Gene model: exons 7, 8, 9 plus the intron 7/8 span
build_gff <- function() {
  rows <- data.frame(
    chr=ZOOM_CHR,
    start=c(EXON7_START, INTRON7_START, EXON8_START, INTRON8_START, EXON9_START),
    end  =c(EXON7_END,   INTRON7_END,   EXON8_END,   INTRON8_END,   EXON9_END),
    type =c("exon","intron","exon","intron","exon"),
    name =c("SMN2_E7","SMN2_intron7","SMN2_E8","SMN2_intron8","SMN2_E9"),
    stringsAsFactors=FALSE
  )
  GRanges(rows$chr, IRanges(rows$start, rows$end), strand="+",
          type=rows$type, name=rows$name)
}
GEs <- build_gff()

m_aso_ctrl <- pooled$ASO_CTRL; m_aso_ctrl <- m_aso_ctrl[m_aso_ctrl$readsN >= 4]
m_scr_ctrl <- pooled$Scramble_CTRL; m_scr_ctrl <- m_scr_ctrl[m_scr_ctrl$readsN >= 4]
m_aso_vpa  <- pooled$ASO_VPA; m_aso_vpa <- m_aso_vpa[m_aso_vpa$readsN >= 4]
m_scr_vpa  <- pooled$Scramble_VPA; m_scr_vpa <- m_scr_vpa[m_scr_vpa$readsN >= 4]

region <- GRanges(ZOOM_CHR, IRanges(PLOT_START, PLOT_END))

plot_col_aso  <- unname(c(COND_COLOURS["ASO_CTRL"], COND_COLOURS["Scramble_CTRL"],
                          "#2D2D2D", "#999999", "#27AE60"))
plot_col_combo <- unname(c(COND_COLOURS["ASO_VPA"], COND_COLOURS["Scramble_CTRL"],
                           "#2D2D2D", "#999999", "#27AE60"))

# overlay enhancer/cCRE annotation as a "DMR-like" marker if any were found
annot_overlay <- if (length(enh_in_window) > 0 || nrow(ccre_enh) > 0) {
  combined <- c(enh_in_window, ccre_gr[grepl("ELS", ccre_sub$type)])
  mcols(combined)$regionType <- "gain"
  combined
} else {
  GRanges()
}

out_path1 <- file.path(OUT_DIR, "SMN2_introns7_8_enhancer_check_ASO_alone.pdf")
pdf(out_path1, width=12, height=6.5)
par(mar=c(5,4,3,1)+0.1, cex=0.9, bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")
plotLocalMethylationProfile(
  methylationData1 = m_aso_ctrl, methylationData2 = m_scr_ctrl,
  region = region,
  DMRs = if (length(annot_overlay)>0) list("Enhancer/cCRE"=annot_overlay) else NULL,
  conditionsNames = c("ASO_CTRL", "Scramble_CTRL"),
  gff = GEs, windowSize = 100, context = "CG", col = plot_col_aso,
  main = "SMN2 introns 7-8: ASO alone vs Scramble CTRL (enhancer check)",
  plotMeanLines = TRUE, plotPoints = TRUE
)
text(par("usr")[1], par("usr")[3]+(par("usr")[4]-par("usr")[3])*0.05,
     labels="Exon 7 (ASO target) | Intron 7 | Exon 8 | Intron 8 | Exon 9",
     cex=0.65, font=3, adj=c(0,0.5))
dev.off()
message("Saved: ", out_path1)

out_path2 <- file.path(OUT_DIR, "SMN2_introns7_8_enhancer_check_combination.pdf")
pdf(out_path2, width=12, height=6.5)
par(mar=c(5,4,3,1)+0.1, cex=0.9, bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")
plotLocalMethylationProfile(
  methylationData1 = m_aso_vpa, methylationData2 = m_scr_ctrl,
  region = region,
  DMRs = if (length(annot_overlay)>0) list("Enhancer/cCRE"=annot_overlay) else NULL,
  conditionsNames = c("ASO_VPA", "Scramble_CTRL"),
  gff = GEs, windowSize = 100, context = "CG", col = plot_col_combo,
  main = "SMN2 introns 7-8: combination vs Scramble CTRL (enhancer check)",
  plotMeanLines = TRUE, plotPoints = TRUE
)
text(par("usr")[1], par("usr")[3]+(par("usr")[4]-par("usr")[3])*0.05,
     labels="Exon 7 (ASO target) | Intron 7 | Exon 8 | Intron 8 | Exon 9",
     cex=0.65, font=3, adj=c(0,0.5))
dev.off()
message("Saved: ", out_path2)

message("\nDone.")
