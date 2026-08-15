.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

# ---- Dedicated, tightly-zoomed locus plot for the SMN2 3' downstream DMR ----
# chr5:70,088,223-70,088,522, ~17.5kb downstream of exon 7, ~10.6kb downstream of the gene's 3 end (gene body ends
# at 70,078,522). DMR called only in the combination contrast at locked
# thresholds (p=6.74e-05, +36.3%). Permutation test on raw proportions:
# deviation=-0.169, p=0.43 -- not significant given low coverage (4-6 CpGs).
# No H9 enhancer or ENCODE cCRE overlap at this position (confirmed earlier).

BY_CHR_UNMASK <- "results/alignments/bs/by_chr"
OUT_DIR       <- "results/smn2_enhancer"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

FLANK    <- 2000
WIN_SIZE <- 100   # smaller window since this is a tight, low-coverage zoom

SMN2_END   <- 70078522
DMR_WINDOW <- list(start=70088223, end=70088522)
PLOT_START <- DMR_WINDOW$start - FLANK
PLOT_END   <- DMR_WINDOW$end + FLANK

COND_COLOURS <- c(
  ASO_CTRL      = "#1B4F8A",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500"
)
GENE_COLOUR <- "#2D2D2D"
TE_COLOUR   <- "#999999"
DMR_COLOUR  <- "#B2182B"

COMPARISONS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond1="ASO_CTRL", cond2="Scramble_CTRL",
       label="ASO effect (CTRL background)"),
  list(name="ASO_VPA_vs_Scramble_CTRL",
       cond1="ASO_VPA", cond2="Scramble_CTRL",
       label="Combined vs CTRL")
)
NEEDED <- unique(unlist(lapply(COMPARISONS, function(x) c(x$cond1, x$cond2))))

# No gene body in this window (it's ~10.6kb downstream of SMN2's 3' end), so the gff
# only marks the SMN2 3' end as a labelled reference point, not a real gene.
build_gff <- function() {
  GRanges(
    seqnames = "chr5",
    ranges   = IRanges(SMN2_END - 200, SMN2_END),
    strand   = "+",
    type     = "gene",
    name     = "SMN2_3prime_end"
  )
}
GEs <- build_gff()

read_unmasked_cpg <- function(condition) {
  message("  loading: ", condition)
  files <- file.path(BY_CHR_UNMASK,
                     sprintf("%s_%d_chr5.CpG_report.txt.gz", condition, 1:3))
  files <- files[file.exists(files)]
  if (length(files) == 0) { message("    no files found"); return(NULL) }
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

message("Loading chr5 methylation for: ", paste(NEEDED, collapse=", "))
pooled <- setNames(lapply(NEEDED, read_unmasked_cpg), NEEDED)

build_dmr_overlay <- function(region_type) {
  GRanges(
    seqnames   = "chr5",
    ranges     = IRanges(DMR_WINDOW$start, DMR_WINDOW$end),
    regionType = region_type
  )
}

for (comp in COMPARISONS) {
  message("\n=== ", comp$name, " ===")
  m1 <- pooled[[comp$cond1]]; m1 <- m1[m1$readsN >= 1]
  m2 <- pooled[[comp$cond2]]; m2 <- m2[m2$readsN >= 1]

  region <- GRanges("chr5", IRanges(PLOT_START, PLOT_END))

  # ASO alone: data showed ASO_CTRL slightly LOWER than Scramble_CTRL at
  # this window (0.939 vs 0.974) -> cond1 lower -> "gain" by convention.
  # Combination: ASO_VPA much lower than Scramble_CTRL (0.596 vs 0.974
  # baseline; vs 0.974 Scramble_CTRL specifically) -> also "gain".
  region_type <- "gain"
  dmrs <- if (comp$name == "ASO_VPA_vs_Scramble_CTRL") build_dmr_overlay(region_type) else GRanges()

  out_path <- file.path(OUT_DIR, paste0("SMN2_3prime_DMR_zoom_", comp$name, ".pdf"))
  plot_col <- unname(c(COND_COLOURS[comp$cond1], COND_COLOURS[comp$cond2],
                        GENE_COLOUR, TE_COLOUR, DMR_COLOUR))
  message("  Using colours: ", paste(plot_col, collapse=", "))

  pdf(out_path, width=10, height=6)
  par(mar=c(5,4,3,1)+0.1, cex=0.9,
      bg="white", col.axis="black", col.lab="black",
      col.main="black", fg="black")
  plotLocalMethylationProfile(
    methylationData1 = m1,
    methylationData2 = m2,
    region           = region,
    DMRs             = if (length(dmrs) > 0) list("DMR"=dmrs) else NULL,
    conditionsNames  = c(comp$cond1, comp$cond2),
    gff              = GEs,
    windowSize       = WIN_SIZE,
    context          = "CG",
    col              = plot_col,
    main             = sprintf("SMN2 3' downstream DMR: %s vs %s (%s)",
                               comp$cond1, comp$cond2, comp$label),
    plotMeanLines    = TRUE,
    plotPoints       = TRUE
  )
  usr <- par("usr")
  text(usr[1], usr[3] + (usr[4]-usr[3])*0.05,
       labels="SMN2's 3 end is ~10.6kb upstream of this window (distal intergenic, nearest-gene call)",
       cex=0.65, font=3, adj=c(0,0.5))
  dev.off()
  message("  Saved: ", basename(out_path))
}

message("\nDone. Outputs in: ", OUT_DIR)
