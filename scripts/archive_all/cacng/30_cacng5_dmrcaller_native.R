.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

# ---- CACNG5 locus methylation profile via DMRcaller native plotting ----
# Uses DMRcaller::plotLocalMethylationProfile -- NOT a custom base-R or
# ggplot2 reimplementation. Modelled directly on 05_smn2_locus_final.R.

BY_CHR_UNMASK <- "results/alignments/bs/by_chr"
OUT_DIR       <- "results/cacng_cluster_dmrcaller"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

FLANK    <- 2000
WIN_SIZE <- 300

# CACNG5 gene + DMR coordinates (GRCh38, chr17)
LOCUS <- list(chr="chr17", start=66893661, end=66960714, strand="+")
DMR_WINDOW <- list(start=66911154, end=66911753)

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

# ---- gene + exon track for DMRcaller plot ----
# CACNG5 exons (GRCh38, chr17) -- matches the pattern in 05_smn2_locus_final.R
EXONS_CACNG5 <- data.frame(
  exon  = c("E1","E2","E3","E4"),
  start = c(66893661, 66910743, 66940481, 66956284),
  end   = c(66893947, 66911770, 66940712, 66960714),
  is_target = c(FALSE, TRUE, FALSE, FALSE)  # E2 overlaps the DMR window
)

build_gff <- function() {
  rows <- list()
  rows[[1]] <- data.frame(
    chr=LOCUS$chr, start=LOCUS$start, end=LOCUS$end,
    strand=LOCUS$strand, type="gene", name="CACNG5", stringsAsFactors=FALSE)
  for (i in seq_len(nrow(EXONS_CACNG5)))
    rows[[length(rows)+1]] <- data.frame(
      chr=LOCUS$chr, start=EXONS_CACNG5$start[i], end=EXONS_CACNG5$end[i],
      strand=LOCUS$strand, type="exon",
      name=sprintf("CACNG5_%s", EXONS_CACNG5$exon[i]), stringsAsFactors=FALSE)
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          strand=df$strand, type=df$type, name=df$name)
}
GEs <- build_gff()

# ---- per-condition CpG reader (pool 3 replicates, unmasked chr17) ----
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

# ---- DMR object for shaded overlay (from confirmed annotated CSV values) ----
# CACNG5 DMR at chr17:66,911,154-66,911,753 in ASO_VPA_vs_Scramble_CTRL,
# regionType "gain" (combination is hypomethylated, confirmed 2026-06-19)
build_dmr_overlay <- function(region_type) {
  GRanges(
    seqnames   = LOCUS$chr,
    ranges     = IRanges(DMR_WINDOW$start, DMR_WINDOW$end),
    regionType = region_type
  )
}

# ---- plot each comparison using DMRcaller's native function ----
for (comp in COMPARISONS) {
  message("\n=== ", comp$name, " ===")
  m1 <- pooled[[comp$cond1]]; m1 <- m1[m1$readsN >= 10]
  m2 <- pooled[[comp$cond2]]; m2 <- m2[m2$readsN >= 10]

  region <- GRanges(seqnames=LOCUS$chr,
                    ranges=IRanges(LOCUS$start-FLANK, LOCUS$end+FLANK))

  # ASO alone: hypermethylated in cond1 (ASO_CTRL) -> regionType "loss"
  # Combination: hypomethylated in cond1 (ASO_VPA) -> regionType "gain"
  region_type <- if (comp$name == "ASO_CTRL_vs_Scramble_CTRL") "loss" else "gain"
  dmrs <- build_dmr_overlay(region_type)

  out_path <- file.path(OUT_DIR, paste0("CACNG5_", comp$name, "_dmrcaller_native.pdf"))
  plot_col <- unname(c(COND_COLOURS[comp$cond1], COND_COLOURS[comp$cond2]))
  message("  Using colours: ", paste(plot_col, collapse=", "))

  pdf(out_path, width=11, height=6)
  par(mar=c(5,4,3,1)+0.1, cex=0.9,
      bg="white", col.axis="black", col.lab="black",
      col.main="black", fg="black")
  plotLocalMethylationProfile(
    methylationData1 = m1,
    methylationData2 = m2,
    region           = region,
    DMRs             = list("DMRs"=dmrs),
    conditionsNames  = c(comp$cond1, comp$cond2),
    gff              = GEs,
    windowSize       = WIN_SIZE,
    context          = "CG",
    col              = plot_col,
    main             = sprintf("CACNG5: %s vs %s (%s)",
                               comp$cond1, comp$cond2, comp$label),
    plotMeanLines    = TRUE,
    plotPoints       = TRUE
  )
  # manual exon labels -- matches pattern in 05_smn2_locus_final.R since
  # plotLocalMethylationProfile's own gff rendering is subtle at this scale
  for (i in seq_len(nrow(EXONS_CACNG5))) {
    mtext(EXONS_CACNG5$exon[i], side=1,
          at=(EXONS_CACNG5$start[i]+EXONS_CACNG5$end[i])/2,
          line=-1.5, cex=0.45,
          col=if(EXONS_CACNG5$is_target[i]) "red" else "black",
          font=if(EXONS_CACNG5$is_target[i]) 2 else 1)
  }
  usr <- par("usr")
  text(usr[1], usr[3] + (usr[4]-usr[3])*0.05,
       labels="CACNG5", cex=0.8, font=2, adj=c(0,0.5))
  dev.off()
  message("  Saved: ", basename(out_path))
}

message("\nDone. Outputs in: ", OUT_DIR)
