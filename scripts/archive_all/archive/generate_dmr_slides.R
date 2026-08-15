#!/usr/bin/env Rscript
# ---
# generate_dmr_slides.R
# Auto-generate a summary slide table from DMR results RDS files.
# Reads the three contrast RDS files and prints a formatted summary table
# that can be pasted into the meeting slides.
#
# Also writes a TSV that can be imported into PowerPoint directly.
#
# Usage: Rscript scripts/generate_dmr_slides.R
# (run from project root after dmrcaller_combine_chr.R completes)
# ---
 
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
 
.libPaths(c("~/R/library", .libPaths()))
 
DMR_DIR <- "results/dmr"
OUT_DIR <- "results/qc/slides"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
 
# Define contrasts and their plain-English labels
CONTRASTS <- list(
  list(
    rds    = "dmr_ASO_VPA_vs_Scramble_CTRL.rds",
    label  = "ASO+VPA vs Scramble_CTRL",
    role   = "Primary (combination vs baseline)"
  ),
  list(
    rds    = "dmr_Scramble_VPA_vs_Scramble_CTRL.rds",
    label  = "VPA alone vs Scramble_CTRL",
    role   = "VPA effect (HDAC inhibitor)"
  ),
  list(
    rds    = "dmr_ASO_CTRL_vs_Scramble_CTRL.rds",
    label  = "ASO alone vs Scramble_CTRL",
    role   = "Negative control (nusinersen only)"
  )
)
 
# Build summary table
rows <- lapply(CONTRASTS, function(ct) {
  rds_path <- file.path(DMR_DIR, ct$rds)
  if (!file.exists(rds_path)) {
    message("Missing: ", rds_path, " — skipping")
    return(NULL)
  }
  dmrs    <- readRDS(rds_path)
  n_total <- length(dmrs)
 
  # regionType convention in DMRcaller v0.25.1:
  # "gain" = proportion1 < proportion2 = hypomethylated in treatment
  # "loss" = proportion1 > proportion2 = hypermethylated in treatment
  n_hypo  <- sum(dmrs$regionType == "gain", na.rm = TRUE)
  n_hyper <- sum(dmrs$regionType == "loss", na.rm = TRUE)
  n_chr13 <- sum(as.character(seqnames(dmrs)) == "chr13", na.rm = TRUE)
 
  # Per-chromosome breakdown — useful for checking chr13 enrichment
  chr_counts <- sort(table(as.character(seqnames(dmrs))), decreasing = TRUE)
  top3_chr   <- paste(names(chr_counts)[1:min(3, length(chr_counts))],
                      chr_counts[1:min(3, length(chr_counts))],
                      sep = ":", collapse = ", ")
 
  data.frame(
    contrast      = ct$label,
    role          = ct$role,
    n_DMRs        = format(n_total, big.mark = ","),
    n_hypo        = format(n_hypo,  big.mark = ","),
    n_hyper       = format(n_hyper, big.mark = ","),
    pct_hypo      = if (n_total > 0) paste0(round(100 * n_hypo / n_total, 1), "%") else "NA",
    n_chr13       = format(n_chr13, big.mark = ","),
    pct_chr13     = if (n_total > 0) paste0(round(100 * n_chr13 / n_total, 1), "%") else "NA",
    top3_chr      = top3_chr,
    stringsAsFactors = FALSE
  )
})
 
rows    <- Filter(Negate(is.null), rows)
summary <- do.call(rbind, rows)
 
# Print formatted table to console
cat("\n========================================\n")
cat("  DMR RESULTS SUMMARY (for slides)\n")
cat("========================================\n\n")
cat(sprintf("%-45s %10s %10s %10s %8s %8s %8s\n",
            "Contrast", "Total DMRs", "Hypo", "Hyper", "% Hypo", "Chr13", "% Chr13"))
cat(strrep("-", 105), "\n")
for (i in seq_len(nrow(summary))) {
  cat(sprintf("%-45s %10s %10s %10s %8s %8s %8s\n",
              summary$contrast[i],
              summary$n_DMRs[i],
              summary$n_hypo[i],
              summary$n_hyper[i],
              summary$pct_hypo[i],
              summary$n_chr13[i],
              summary$pct_chr13[i]))
}
cat("\nTop 3 chromosomes by DMR count:\n")
for (i in seq_len(nrow(summary))) {
  cat(sprintf("  %-35s %s\n", summary$contrast[i], summary$top3_chr[i]))
}
 
# Save TSV for PowerPoint import
# This TSV can be opened in Excel and pasted into PowerPoint as a table.
tsv_path <- file.path(OUT_DIR, "dmr_results_for_slides.tsv")
write.table(summary, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("\nSaved slide table: ", tsv_path)
 
# Per-chromosome breakdown
# Useful for checking the chr13 hotspot and identifying other enriched chromosomes
message("\n=== Per-chromosome DMR counts ===")
chr_summary_list <- lapply(CONTRASTS, function(ct) {
  rds_path <- file.path(DMR_DIR, ct$rds)
  if (!file.exists(rds_path)) return(NULL)
  dmrs    <- readRDS(rds_path)
  chr_tab <- as.data.frame(table(chr = as.character(seqnames(dmrs))))
  chr_tab$contrast <- ct$label
  chr_tab <- chr_tab[order(-chr_tab$Freq), ]
  chr_tab
})
chr_summary_list <- Filter(Negate(is.null), chr_summary_list)
chr_summary      <- do.call(rbind, chr_summary_list)
 
chr_tsv_path <- file.path(OUT_DIR, "dmr_per_chromosome_for_slides.tsv")
write.table(chr_summary, chr_tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved per-chromosome breakdown: ", chr_tsv_path)
 
message("\nDone.")
