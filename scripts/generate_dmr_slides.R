#!/usr/bin/env Rscript
# =============================================================================
# generate_dmr_slides.R
# Pull numbers from DMR result RDS files and format them for slides.
# Writes a formatted console table + TSV you can paste into PowerPoint.
#
# Usage: Rscript scripts/generate_dmr_slides.R
# Run from project root after dmrcaller_combine_chr.R completes.
# =============================================================================
 
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
 
.libPaths(c("~/R/library", .libPaths()))
 
DMR_DIR <- "results/dmr"
OUT_DIR <- "results/qc/slides"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
 
CONTRASTS <- list(
  list(rds = "dmr_ASO_VPA_vs_Scramble_CTRL.rds",
       label = "ASO+VPA vs Scramble_CTRL",
       role  = "Primary (combination vs baseline)"),
  list(rds = "dmr_Scramble_VPA_vs_Scramble_CTRL.rds",
       label = "VPA alone vs Scramble_CTRL",
       role  = "VPA effect (HDAC inhibitor)"),
  list(rds = "dmr_ASO_CTRL_vs_Scramble_CTRL.rds",
       label = "ASO alone vs Scramble_CTRL",
       role  = "Negative control (nusinersen only)")
)
 
rows <- lapply(CONTRASTS, function(ct) {
  rds_path <- file.path(DMR_DIR, ct$rds)
  if (!file.exists(rds_path)) { message("Missing: ", rds_path); return(NULL) }
  dmrs    <- readRDS(rds_path)
  n_total <- length(dmrs)
 
  # "gain" = hypomethylated in treatment, "loss" = hypermethylated
  n_hypo  <- sum(dmrs$regionType == "gain", na.rm = TRUE)
  n_hyper <- sum(dmrs$regionType == "loss", na.rm = TRUE)
  n_chr13 <- sum(as.character(seqnames(dmrs)) == "chr13", na.rm = TRUE)
 
  chr_counts <- sort(table(as.character(seqnames(dmrs))), decreasing = TRUE)
  top3 <- paste(names(chr_counts)[1:min(3, length(chr_counts))],
                chr_counts[1:min(3, length(chr_counts))],
                sep = ":", collapse = ", ")
 
  data.frame(
    contrast  = ct$label,
    role      = ct$role,
    n_DMRs    = format(n_total, big.mark = ","),
    n_hypo    = format(n_hypo,  big.mark = ","),
    n_hyper   = format(n_hyper, big.mark = ","),
    pct_hypo  = if (n_total > 0) paste0(round(100 * n_hypo  / n_total, 1), "%") else "NA",
    n_chr13   = format(n_chr13, big.mark = ","),
    pct_chr13 = if (n_total > 0) paste0(round(100 * n_chr13 / n_total, 1), "%") else "NA",
    top3_chr  = top3,
    stringsAsFactors = FALSE
  )
})
 
rows    <- Filter(Negate(is.null), rows)
summary <- do.call(rbind, rows)
 
cat("\n===================================================\n")
cat("  DMR RESULTS SUMMARY\n")
cat("===================================================\n\n")
cat(sprintf("%-40s %10s %10s %10s %8s %8s %8s\n",
            "Contrast", "Total", "Hypo", "Hyper", "% Hypo", "Chr13", "% Chr13"))
cat(strrep("-", 100), "\n")
for (i in seq_len(nrow(summary))) {
  cat(sprintf("%-40s %10s %10s %10s %8s %8s %8s\n",
              summary$contrast[i], summary$n_DMRs[i],
              summary$n_hypo[i],   summary$n_hyper[i],
              summary$pct_hypo[i], summary$n_chr13[i],
              summary$pct_chr13[i]))
}
 
cat("\nTop 3 chromosomes by DMR count:\n")
for (i in seq_len(nrow(summary))) {
  cat(sprintf("  %-35s %s\n", summary$contrast[i], summary$top3_chr[i]))
}
 
write.table(summary,
            file.path(OUT_DIR, "dmr_results_for_slides.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("\nSaved: ", file.path(OUT_DIR, "dmr_results_for_slides.tsv"))
 
# Per-chromosome breakdown
chr_rows <- lapply(CONTRASTS, function(ct) {
  rds_path <- file.path(DMR_DIR, ct$rds)
  if (!file.exists(rds_path)) return(NULL)
  dmrs    <- readRDS(rds_path)
  chr_tab <- as.data.frame(table(chr = as.character(seqnames(dmrs))))
  chr_tab$contrast <- ct$label
  chr_tab[order(-chr_tab$Freq), ]
})
chr_rows <- Filter(Negate(is.null), chr_rows)
write.table(do.call(rbind, chr_rows),
            file.path(OUT_DIR, "dmr_per_chromosome.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: ", file.path(OUT_DIR, "dmr_per_chromosome.tsv"))
message("\nDone.")
