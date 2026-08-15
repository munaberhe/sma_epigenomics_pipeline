.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

# ---- Distal intergenic DMR downstream of CACNG5 (chr17), DMRcaller native plot ----
# v4 fixes:
#   1. col vector now has >= 4 elements (cond1, cond2, gene, TE) as required by
#      DMRcaller's .isColor(col, minLength = 4 + numberOfDMRs) check -- a 2-color
#      vector silently failed validation and fell back to DMRcaller's own
#      default palette, which is why colors were wrong in v1-v3.
#   2. Real CACNG5/CACNG4/PRKCA coordinates from TxDb.Hsapiens.UCSC.hg38.knownGene,
#      confirmed 2026-06-20: CACNG5 66,835,117-66,894,751; CACNG4 66,964,707-67,033,398;
#      PRKCA 66,302,613-66,810,743.
#   3. DMR at chr17:66,911,154-66,911,753 is Distal Intergenic (ChIPseeker),
#      distanceToTSS=33,881 from CACNG5 -- nearest-gene annotation, NOT intragenic.
#      Plot window widened to show CACNG5's real gene body alongside the DMR
#      so the intergenic position is visually obvious, not implied otherwise.

BY_CHR_UNMASK <- "results/alignments/bs/by_chr"
OUT_DIR       <- "results/cacng_cluster_dmrcaller"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

FLANK    <- 3000
WIN_SIZE <- 300

# Real coordinates, confirmed via TxDb.Hsapiens.UCSC.hg38.knownGene
CACNG5 <- list(chr="chr17", start=66835117, end=66894751, strand="+")
DMR_WINDOW <- list(start=66911154, end=66911753)

# Window spans CACNG5's real gene body through the DMR position
PLOT_START <- CACNG5$start - FLANK
PLOT_END   <- DMR_WINDOW$end + FLANK

COMPARISONS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond1="ASO_CTRL", cond2="Scramble_CTRL",
       label="ASO effect (CTRL background)"),
  list(name="ASO_VPA_vs_Scramble_CTRL",
       cond1="ASO_VPA", cond2="Scramble_CTRL",
       label="Combined vs CTRL")
)

NEEDED <- unique(unlist(lapply(COMPARISONS, function(x) c(x$cond1, x$cond2))))

COND_COLOURS <- c(
  ASO_CTRL      = "#1B4F8A",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500"
)
GENE_COLOUR <- "#2D2D2D"
TE_COLOUR   <- "#999999"
DMR_COLOUR  <- "#B2182B"

# ---- gene track: real CACNG5 exons from TxDb ----
build_gff <- function() {
  rows <- list()
  rows[[1]] <- data.frame(
    chr=CACNG5$chr, start=CACNG5$start, end=CACNG5$end,
    strand=CACNG5$strand, type="gene", name="CACNG5", stringsAsFactors=FALSE)
  # Real exon coordinates from TxDb (chr17, gene CACNG5, ENTREZID 27091)
  exon_starts <- c(66835117,66877230,66877269,66877273,66878972,
                   66880557,66884516,66884983,66884983,66893412)
  exon_ends   <- c(66835250,66877528,66877528,66877528,66879058,
                   66880697,66884661,66885278,66894751,66894742)
  for (i in seq_along(exon_starts))
    rows[[length(rows)+1]] <- data.frame(
      chr=CACNG5$chr, start=exon_starts[i], end=exon_ends[i],
      strand=CACNG5$strand, type="exon",
      name=sprintf("CACNG5_E%d", i), stringsAsFactors=FALSE)
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          strand=df$strand, type=df$type, name=df$name)
}
GEs <- build_gff()

read_unmasked_cpg <- function(condition) {
  message("  loading: ", condition)
  files <- file.path(BY_CHR_UNMASK,
                     sprintf("%s_%d_chr17.CpG_report.txt.gz", condition, 1:3))
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

message("Loading chr17 methylation for: ", paste(NEEDED, collapse=", "))
pooled <- setNames(lapply(NEEDED, read_unmasked_cpg), NEEDED)

build_dmr_overlay <- function(region_type) {
  GRanges(
    seqnames   = CACNG5$chr,
    ranges     = IRanges(DMR_WINDOW$start, DMR_WINDOW$end),
    regionType = region_type
  )
}

for (comp in COMPARISONS) {
  message("\n=== ", comp$name, " ===")
  m1 <- pooled[[comp$cond1]]; m1 <- m1[m1$readsN >= 10]
  m2 <- pooled[[comp$cond2]]; m2 <- m2[m2$readsN >= 10]

  region <- GRanges(seqnames=CACNG5$chr, ranges=IRanges(PLOT_START, PLOT_END))

  region_type <- if (comp$name == "ASO_CTRL_vs_Scramble_CTRL") "loss" else "gain"
  dmrs <- build_dmr_overlay(region_type)

  out_path <- file.path(OUT_DIR, paste0("CACNG5_region_", comp$name, "_dmrcaller_v4.pdf"))

  # Color vector: cond1, cond2, gene, TE, DMR (>= 4 + numberOfDMRs required)
  plot_col <- unname(c(COND_COLOURS[comp$cond1], COND_COLOURS[comp$cond2],
                        GENE_COLOUR, TE_COLOUR, DMR_COLOUR))
  message("  Using colours: ", paste(plot_col, collapse=", "))

  pdf(out_path, width=12, height=6.5)
  par(mar=c(5,4,3,1)+0.1, cex=0.9,
      bg="white", col.axis="black", col.lab="black",
      col.main="black", fg="black")
  plotLocalMethylationProfile(
    methylationData1 = m1,
    methylationData2 = m2,
    region           = region,
    DMRs             = list("DMR"=dmrs),
    conditionsNames  = c(comp$cond1, comp$cond2),
    gff              = GEs,
    windowSize       = WIN_SIZE,
    context          = "CG",
    col              = plot_col,
    main             = sprintf("Distal intergenic DMR near CACNG5: %s vs %s (%s)",
                               comp$cond1, comp$cond2, comp$label),
    plotMeanLines    = TRUE,
    plotPoints       = TRUE
  )
  usr <- par("usr")
  text(usr[1], usr[3] + (usr[4]-usr[3])*0.05,
       labels="CACNG5 (gene body)  |  DMR is 33.9 kb downstream (distal intergenic)",
       cex=0.65, font=3, adj=c(0,0.5))
  dev.off()
  message("  Saved: ", basename(out_path))
}

message("\nDone. Outputs in: ", OUT_DIR)
