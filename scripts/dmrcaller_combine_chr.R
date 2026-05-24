#!/usr/bin/env Rscript
# =============================================================================
# dmrcaller_combine_chr.R
# Combine per-chromosome DMR results into genome-wide GRanges objects.
# SMA Epigenomics Pipeline — Muna Berhe · bt25018 · QMUL
#
# Run after all 72 per-chromosome jobs (dmrcaller_by_chr.R) complete.
# Takes 24 per-chromosome RDS files per contrast and concatenates them into
# one genome-wide GRanges, then writes RDS + BED outputs.
#
# Can also be run on partial results (e.g. if some chr jobs are still running)
# to get a preview — missing chromosomes are skipped with a warning.
#
# Usage: Rscript scripts/dmrcaller_combine_chr.R
# =============================================================================
 
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
 
.libPaths(c("~/R/library", .libPaths()))
 
CHROMS <- paste0("chr", c(1:22, "X", "Y"))
 
CONTRASTS <- c(
  "ASO_VPA_vs_Scramble_CTRL",       # combination vs baseline — primary result
  "Scramble_VPA_vs_Scramble_CTRL",  # VPA alone — isolates HDAC inhibitor effect
  "ASO_CTRL_vs_Scramble_CTRL"       # ASO alone — negative control
)
 
OUT_DIR <- "results/dmr"
 
for (contrast in CONTRASTS) {
  message("\nCombining: ", contrast)
  chr_dir <- file.path(OUT_DIR, "by_chr", contrast)
 
  # Load per-chromosome RDS files. NULL returned for missing chromosomes.
  # Filter these out below so partial results still produce a usable output.
  dmr_list <- lapply(CHROMS, function(chr) {
    f <- file.path(chr_dir, paste0("dmr_", contrast, "_", chr, ".rds"))
    if (!file.exists(f)) {
      message("  Missing: ", chr, " — job may still be running")
      return(NULL)
    }
    readRDS(f)
  })
 
  dmr_list <- Filter(Negate(is.null), dmr_list)           # remove missing chromosomes
  dmr_list <- Filter(function(x) length(x) > 0, dmr_list) # remove chromosomes with 0 DMRs
 
  if (length(dmr_list) == 0) { message("  No results"); next }
 
  # do.call(c, ...) concatenates GRanges objects — equivalent to rbind for data frames.
  # Produces one genome-wide GRanges with DMRs from all chromosomes.
  all_dmrs <- do.call(c, dmr_list)
  n_total  <- length(all_dmrs)
 
  # regionType convention in DMRcaller v0.25.1:
  # "gain" = proportion1 < proportion2 = hypomethylated in treatment (cond_a)
  # "loss" = proportion1 > proportion2 = hypermethylated in treatment (cond_a)
  # Confirmed empirically: chr13 DMRs have proportion1 (ASO_VPA) ~0.27-0.55
  # vs proportion2 (Scramble_CTRL) ~0.48-0.78 but regionType = "gain".
  n_hypo  <- sum(all_dmrs$regionType == "gain", na.rm = TRUE)
  n_hyper <- sum(all_dmrs$regionType == "loss", na.rm = TRUE)
 
  message("  Total DMRs: ", n_total, " (", n_hyper, " hyper / ", n_hypo, " hypo)")
  message("  Chr13: ", sum(as.character(seqnames(all_dmrs)) == "chr13"))
  message("  Chromosomes with results: ", length(dmr_list), "/", length(CHROMS))
 
  # Save RDS — primary input for dmr_annotate.R (ChIPseeker + clusterProfiler)
  saveRDS(all_dmrs, file.path(OUT_DIR, paste0("dmr_", contrast, ".rds")))
 
  # BED file for IGV / genome browser
  # GRanges is 1-based; BED is 0-based half-open, so start = GRanges_start - 1
  # score = -log10(pValue): larger = more significant
  # methDiff: +1 = hypermethylated in treatment, -1 = hypomethylated
  bed_df <- data.frame(
    chr      = as.character(seqnames(all_dmrs)),
    start    = start(all_dmrs) - 1,
    end      = end(all_dmrs),
    name     = contrast,
    score    = round(-log10(all_dmrs$pValue + 1e-300), 2),  # 1e-300 prevents log10(0)
    strand   = ".",
    methDiff = round(as.integer(all_dmrs$regionType == "loss") -
                     as.integer(all_dmrs$regionType == "gain"), 4),
    type     = all_dmrs$regionType,
    nCpG     = all_dmrs$cytosinesCount,
    pValue   = all_dmrs$pValue
  )
  write.table(bed_df, file.path(OUT_DIR, paste0("dmr_", contrast, ".bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
 
  message("  Saved RDS and BED to: ", OUT_DIR)
}
 
message("\nDone.")
