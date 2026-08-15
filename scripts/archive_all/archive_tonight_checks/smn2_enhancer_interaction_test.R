.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

BY_CHR_UNMASK <- "results/alignments/bs/by_chr"
WINDOW <- list(chr="chr5", start=70088223, end=70088522)

CONDITIONS <- c("ASO_CTRL", "Scramble_CTRL", "ASO_VPA", "Scramble_VPA")

read_unmasked_cpg <- function(condition, chr) {
  files <- file.path(BY_CHR_UNMASK,
                     sprintf("%s_%d_%s.CpG_report.txt.gz", condition, 1:3, chr))
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

message("Loading chr9 methylation for all four conditions...")
region <- GRanges(WINDOW$chr, IRanges(WINDOW$start, WINDOW$end))

results <- list()
for (cond in CONDITIONS) {
  message("  ", cond)
  pooled <- read_unmasked_cpg(cond, WINDOW$chr)
  if (is.null(pooled)) {
    results[[cond]] <- NA
    next
  }
  hits <- subsetByOverlaps(pooled, region)
  hits <- hits[hits$readsN >= 4]
  total_M <- sum(hits$readsM)
  total_N <- sum(hits$readsN)
  prop <- if (total_N > 0) total_M / total_N else NA
  results[[cond]] <- prop
  message(sprintf("    n CpGs=%d  total reads=%d  proportion=%.4f", length(hits), total_N, prop))
}

cat("\n=== SIGMAR1 window methylation, all four conditions ===\n")
cat(sprintf("ASO_CTRL:      %.4f\n", results$ASO_CTRL))
cat(sprintf("Scramble_CTRL: %.4f\n", results$Scramble_CTRL))
cat(sprintf("ASO_VPA:       %.4f\n", results$ASO_VPA))
cat(sprintf("Scramble_VPA:  %.4f\n", results$Scramble_VPA))

cat("\n=== Additive prediction check ===\n")
aso_effect <- results$ASO_CTRL - results$Scramble_CTRL
vpa_effect <- results$Scramble_VPA - results$Scramble_CTRL
predicted_additive <- results$Scramble_CTRL + aso_effect + vpa_effect
actual_combo <- results$ASO_VPA

cat(sprintf("ASO-alone effect (ASO_CTRL - Scramble_CTRL):   %+.4f\n", aso_effect))
cat(sprintf("VPA-alone effect (Scramble_VPA - Scramble_CTRL): %+.4f\n", vpa_effect))
cat(sprintf("Predicted additive combo value:                  %.4f\n", predicted_additive))
cat(sprintf("Actual observed combo value:                     %.4f\n", actual_combo))
cat(sprintf("Deviation from additive prediction:              %+.4f\n", actual_combo - predicted_additive))
