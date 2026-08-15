# Adapted from smn_locus_dmrcaller_comparisons.R
# Uses SMN1-masked chr5 CX reports from results/alignments_smn1_masked/chr5_cx/
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

.libPaths(c("~/R/library", .libPaths()))
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

# ── Config ────────────────────────────────────────────────────────────────────
CHR5_DIR <- "results/alignments_smn1_masked/chr5_cx"
OUT_DIR  <- "results/smn2_masked_profile"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

FLANK    <- 5000
WIN_SIZE <- 300

LOCI <- list(
  SMN1 = list(chr="chr5", start=70924941, end=70953015, strand="+"),
  SMN2 = list(chr="chr5", start=70049638, end=70078522, strand="+")
)

# Exon coordinates using Alberto's convention
# GTF exon 8 = Alberto E7 = ASO target
EXONS <- list(
  SMN1 = data.frame(
    exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start = c(70925030,70938807,70941357,70942326,70942686,
              70944627,70946033,70951913,70952411),
    end   = c(70925158,70938878,70941476,70942526,70942838,
              70944722,70946143,70951966,70952984),
    is_target = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)
  ),
  SMN2 = data.frame(
    exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start = c(70049638,70063415,70065965,70066934,70067294,
              70069235,70070641,70076521,70077019),
    end   = c(70049766,70063486,70066084,70067134,70067446,
              70069330,70070751,70076574,70077592),
    is_target = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)
  )
)

COMPARISONS <- list(
  list(name="ASO_vs_Scramble_CTRL",
       cond1="ASO_CTRL",  cond2="Scramble_CTRL",
       label="ASO effect (CTRL background)"),
  list(name="VPA_vs_CTRL_ASO",
       cond1="ASO_CTRL",  cond2="ASO_VPA",
       label="VPA effect (ASO background)"),
  list(name="ASO_vs_Scramble_VPA",
       cond1="ASO_VPA",   cond2="Scramble_VPA",
       label="ASO effect (VPA background)")
)

NEEDED_CONDS <- unique(unlist(lapply(COMPARISONS, function(x) c(x$cond1, x$cond2))))

# ── Reader — reads pooled uncompressed chr5 CX reports ───────────────────────
read_condition_masked <- function(condition) {
  reps <- c(1,2,3)
  grs  <- lapply(reps, function(r) {
    sample <- paste0(condition, "_", r)
    path   <- file.path(CHR5_DIR, paste0(sample, "_chr5.CX_report.txt"))
    if (!file.exists(path)) stop("Missing: ", path)
    message("    Reading: ", sample)
    d <- read.table(path, header=FALSE, sep="\t",
                    col.names=c("chr","pos","strand","countM","countU","context","tri"),
                    colClasses=c("character","integer","character",
                                 "integer","integer","character","character"))
    d <- d[d$context == "CG", ]
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos, d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

message("Reading and pooling masked chr5 CX reports...")
pooled_data <- lapply(NEEDED_CONDS, function(c) {
  message("  ", c)
  read_condition_masked(c)
})
names(pooled_data) <- NEEDED_CONDS

# ── Build GFF with Alberto's exon labels ─────────────────────────────────────
build_gff <- function() {
  rows <- list()
  for (locus_name in names(LOCI)) {
    locus <- LOCI[[locus_name]]
    rows[[length(rows)+1]] <- data.frame(
      chr=locus$chr, start=locus$start, end=locus$end,
      strand=locus$strand, type="gene", name=locus_name,
      stringsAsFactors=FALSE)
    ex <- EXONS[[locus_name]]
    for (i in seq_len(nrow(ex))) {
      rows[[length(rows)+1]] <- data.frame(
        chr=locus$chr, start=ex$start[i], end=ex$end[i],
        strand=locus$strand, type="exon",
        name=sprintf("%s_%s", locus_name, ex$exon[i]),
        stringsAsFactors=FALSE)
    }
  }
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          strand=df$strand, type=df$type, name=df$name)
}

GEs <- build_gff()
message("Built GFF: ", length(GEs), " features")

# ── Plot function — identical to original ────────────────────────────────────
plot_one <- function(comp, locus_name) {
  locus  <- LOCI[[locus_name]]
  region <- GRanges(seqnames=locus$chr,
                    ranges=IRanges(locus$start - FLANK, locus$end + FLANK))

  plotLocalMethylationProfile(
    methylationData1 = pooled_data[[comp$cond1]],
    methylationData2 = pooled_data[[comp$cond2]],
    region           = region,
    DMRs             = NULL,
    conditionsNames  = c(comp$cond1, comp$cond2),
    gff              = GEs,
    windowSize       = WIN_SIZE,
    context          = "CG",
    main             = sprintf("%s ... %s vs %s (%s)",
                               locus_name, comp$cond1, comp$cond2, comp$label),
    plotMeanLines    = TRUE,
    plotPoints       = TRUE
  )

  # Exon labels using Alberto's convention
  ex <- EXONS[[locus_name]]
  for (i in seq_len(nrow(ex))) {
    x_mid  <- (ex$start[i] + ex$end[i]) / 2
    colour <- if (ex$is_target[i]) "red" else "black"
    weight <- if (ex$is_target[i]) 2 else 1
    mtext(ex$exon[i], side=1, at=x_mid, line=-1.5,
          cex=0.45, col=colour, font=weight)
  }
}

# ── Individual PDFs ───────────────────────────────────────────────────────────
for (comp in COMPARISONS) {
  pdf_path <- file.path(OUT_DIR,
    sprintf("SMN2_masked_dmrcaller_%s.pdf", comp$name))
  message("\nPlotting: ", comp$name)
  pdf(pdf_path, width=11, height=8.5)
  par(mfrow=c(2,1), mar=c(5,4,3,1)+0.1, cex=0.9)
  plot_one(comp, "SMN1")
  plot_one(comp, "SMN2")
  dev.off()
  message("  Saved: ", pdf_path)
}

# ── Combined PDF — all 3 comparisons ─────────────────────────────────────────
combined_pdf <- file.path(OUT_DIR, "SMN2_masked_all_comparisons.pdf")
message("\nPlotting combined PDF...")
pdf(combined_pdf, width=16, height=14)
par(mfrow=c(3,2), mar=c(5,4,3,1)+0.1, cex=0.7)
for (comp in COMPARISONS) {
  plot_one(comp, "SMN1")
  plot_one(comp, "SMN2")
}
dev.off()
message("Saved: ", combined_pdf)

# ── Weighted mean summary ─────────────────────────────────────────────────────
message("\n=== Methylation summary ===")
summarise <- function(condition_name, locus_name) {
  locus <- LOCI[[locus_name]]
  gr    <- pooled_data[[condition_name]]
  sel   <- as.character(seqnames(gr)) == locus$chr &
           start(gr) >= locus$start & start(gr) <= locus$end &
           gr$readsN > 0
  g <- gr[sel]
  if (length(g)==0) return(data.frame(condition=condition_name,
    locus=locus_name, n_cpg=0L, weighted_mean=NA_real_))
  data.frame(condition=condition_name, locus=locus_name,
             n_cpg=length(g),
             weighted_mean=round(sum(g$readsM)/sum(g$readsN),4))
}

summary_rows <- list()
for (cond in NEEDED_CONDS)
  for (locus_name in names(LOCI))
    summary_rows[[length(summary_rows)+1]] <- summarise(cond, locus_name)

summary_df <- do.call(rbind, summary_rows)
tsv_path <- file.path(OUT_DIR, "SMN_masked_weighted_mean_methylation_v3.tsv")
write.table(summary_df, tsv_path, sep="\t", quote=FALSE, row.names=FALSE)
message("\nSummary:")
print(summary_df, row.names=FALSE)
message("\nDone. Outputs in ", OUT_DIR)
