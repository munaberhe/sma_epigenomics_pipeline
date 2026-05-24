#!/usr/bin/env Rscript
# ---
# dmrcaller_combine_chr.R
# Combine per-chromosome DMR results into genome-wide GRanges objects
# SMA Epigenomics Pipeline — Muna Berhe · bt25018 · QMUL
#
# Run this after all 72 per-chromosome jobs (dmrcaller_by_chr.R) complete.
# Takes the 24 per-chromosome RDS files per contrast and concatenates them
# into a single genome-wide GRanges, then writes RDS + BED outputs.
#
# Usage: Rscript scripts/dmrcaller_combine_chr.R
# (run from project root after submit_dmr_by_chr.sh jobs finish)
# ---
 
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
 
.libPaths(c("~/R/library", .libPaths()))
 
# All 24 standard chromosomes — same list used in dmrcaller_by_chr.R and
# 04_split_by_chr.sh. Must match exactly or files won't be found.
# Unplaced contigs and patches are excluded — low coverage, inflate counts.
CHROMS <- paste0("chr", c(1:22, "X", "Y"))
 
CONTRASTS <- c(
  "ASO_VPA_vs_Scramble_CTRL",        # combination vs baseline — primary result
  "Scramble_VPA_vs_Scramble_CTRL",   # VPA alone — isolates HDAC inhibitor effect
  "ASO_CTRL_vs_Scramble_CTRL"        # ASO alone — negative control
)
 
OUT_DIR <- "results/dmr"
 
for (contrast in CONTRASTS) {
  message("\nCombining: ", contrast)
  chr_dir <- file.path(OUT_DIR, "by_chr", contrast)
 
  # Load each per-chromosome RDS. NULL is returned for missing chromosomes
  # (e.g. if a job failed or is still running) — these are filtered out below.
  # This means you can run combine_chr.R on partial results to get a preview,
  # then rerun after all jobs finish to get the complete genome-wide result.
  dmr_list <- lapply(CHROMS, function(chr) {
    f <- file.path(chr_dir, paste0("dmr_", contrast, "_", chr, ".rds"))
    if (!file.exists(f)) {
      message("  Missing: ", chr, " — skipping (job may still be running)")
      return(NULL)
    }
    readRDS(f)
  })
 
  # Remove NULLs (missing chromosomes) and empty GRanges (chromosomes with 0 DMRs)
  dmr_list <- Filter(Negate(is.null), dmr_list)
  dmr_list <- Filter(function(x) length(x) > 0, dmr_list)
 
  if (length(dmr_list) == 0) {
    message("  No results found for any chromosome — skipping")
    next
  }
 
  # do.call(c, ...) concatenates GRanges objects — equivalent to rbind for data frames.
  # This produces one genome-wide GRanges with all DMRs across all chromosomes.
  # Chromosome order follows CHROMS (chr1, chr2, ..., chr22, X, Y).
  all_dmrs <- do.call(c, dmr_list)
  n_total  <- length(all_dmrs)
 
  # Directionality note
  # DMRcaller's regionType labelling is counterintuitive:
  #   "gain" = methylationData2 > methylationData1 = treatment is HYPOMETHYLATED
  #   "loss" = methylationData2 < methylationData1 = treatment is HYPERMETHYLATED
  # This was confirmed by checking proportion1 vs proportion2 in chr13 DMRs:
  #   proportion1 (ASO_VPA) ~0.27-0.55 < proportion2 (Scramble_CTRL) ~0.48-0.78
  #   regionType = "gain" → treatment has LESS methylation = hypomethylated
  # References: DMRcaller v0.25.1 source; confirmed May 2026
  n_hypo  <- sum(all_dmrs$regionType == "gain", na.rm = TRUE)
  n_hyper <- sum(all_dmrs$regionType == "loss", na.rm = TRUE)
 
  message("  Total DMRs: ", n_total,
          " (", n_hyper, " hyper / ", n_hypo, " hypo)")
  message("  Chr13 DMRs: ",
          sum(as.character(seqnames(all_dmrs)) == "chr13"))
  message("  Chromosomes with DMRs: ", length(dmr_list), "/", length(CHROMS))
 
  # Save genome-wide RDS — primary output for dmr_annotate.R (ChIPseeker + GO/KEGG)
  saveRDS(all_dmrs, file.path(OUT_DIR, paste0("dmr_", contrast, ".rds")))
 
  # BED file for IGV / genome browser
  # BED format is 0-based half-open: start = GRanges_start - 1, end = GRanges_end
  # GRanges uses 1-based coordinates (Bioconductor standard) so we subtract 1
  # from start to convert to BED's 0-based system.
  # score = -log10(pValue): larger score = more significant DMR
  # methDiff: +1 = hypermethylated in treatment, -1 = hypomethylated
  # Note: methDiff sign matches the biological direction, not the regionType label
  bed_df <- data.frame(
    chr      = as.character(seqnames(all_dmrs)),
    start    = start(all_dmrs) - 1,   # 1-based GRanges → 0-based BED
    end      = end(all_dmrs),
    name     = contrast,
    score    = round(-log10(all_dmrs$pValue + 1e-300), 2),  # 1e-300 prevents log10(0)
    strand   = ".",
    # methDiff: "gain" (hypo in treatment) gets -1, "loss" (hyper in treatment) gets +1
    methDiff = round(as.integer(all_dmrs$regionType == "loss") -
                     as.integer(all_dmrs$regionType == "gain"), 4),
    type     = all_dmrs$regionType,
    nCpG     = all_dmrs$cytosinesCount,
    pValue   = all_dmrs$pValue
  )
  write.table(bed_df,
              file.path(OUT_DIR, paste0("dmr_", contrast, ".bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
 
  message("  Saved RDS and BED to: ", OUT_DIR)
}
 
message("\nDone. Genome-wide DMR results in: ", OUT_DIR)
