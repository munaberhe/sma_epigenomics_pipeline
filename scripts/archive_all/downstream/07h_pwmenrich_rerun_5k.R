#!/usr/bin/env Rscript
# 07h_pwmenrich_rerun_5k.R
# Re-runs PWMEnrich for all four 5k contrasts and saves:
#   - pwmenrich_results_<contrast>_5k.rds  (MotifEnrichmentResults -- for groupReport)
#   - pwmenrich_report_<contrast>_5k_v2.rds (filtered MotifEnrichmentReport)
#   - pwmenrich_report_<contrast>_5k_v2.pdf (native plot() table, top 20)
#   - pwmenrich_top10_<contrast>_5k_v2.pdf  (native plot() table, top 10)
# Modelled exactly on tf_pwm_ASO_VPA.R
# Muna Berhe -- bt25018 -- QMUL MSc Bioinformatics

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(PWMEnrich)
  library(PWMEnrich.Hsapiens.background)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

OUT_DIR   <- "results/tf_motif"
KEEP_CHRS <- paste0("chr", c(1:22, "X"))
BLACKLIST <- "CENPB|ZNF274|ZNF93|UW[.]Motif"

message("Loading PWM background...")
data(PWMLogn.hg19.MotifDb.Hsap)
bg <- PWMLogn.hg19.MotifDb.Hsap

CONTRASTS <- list(
  list(
    name    = "ASO_VPA_vs_Scramble_CTRL",
    rds     = "results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds",
    label   = "ASO+VPA vs Scramble_CTRL"
  ),
  list(
    name    = "Scramble_VPA_vs_Scramble_CTRL",
    rds     = "results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds",
    label   = "VPA vs Scramble_CTRL"
  ),
  list(
    name    = "ASO_CTRL_vs_Scramble_CTRL",
    rds     = "results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds",
    label   = "ASO_CTRL vs Scramble_CTRL"
  ),
  list(
    name    = "ASO_VPA_vs_ASO_CTRL",
    rds     = "results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds",
    label   = "ASO_VPA vs ASO_CTRL (synergy)"
  )
)

for (cfg in CONTRASTS) {
  message("\n========== ", cfg$label, " ==========")

  if (!file.exists(cfg$rds)) {
    message("  SKIP -- RDS not found: ", cfg$rds)
    next
  }

  # Load and filter DMRs
  gr <- readRDS(cfg$rds)
  gr <- gr[gr$context == "CG" & as.character(seqnames(gr)) %in% KEEP_CHRS]
  message("  DMRs after CG+chr filter: ", length(gr))

  # Sample to 5k if needed
  if (length(gr) > 5000) {
    set.seed(42)
    gr <- gr[sample(length(gr), 5000)]
    message("  Sampled to 5000")
  }
  message("  Final n DMRs: ", length(gr))

  # Extract sequences
  message("  Extracting sequences from hg38...")
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)

  # Run motif enrichment
  message("  Running motifEnrichment (this takes ~10-20 min)...")
  res <- motifEnrichment(seqs, bg, verbose=FALSE)

  # Save raw results RDS -- this is what groupReport() needs
  res_out <- file.path(OUT_DIR,
    paste0("pwmenrich_results_", cfg$name, "_5k.rds"))
  saveRDS(res, res_out)
  message("  Saved results RDS: ", basename(res_out))

  # Generate report
  rep <- groupReport(res)

  # Save report RDS
  rep_out <- file.path(OUT_DIR,
    paste0("pwmenrich_report_", cfg$name, "_5k_v2.rds"))
  saveRDS(rep, rep_out)
  message("  Saved report RDS: ", basename(rep_out))

  # Filter blacklist and deduplicate
  rep_df    <- as.data.frame(rep)
  keep      <- !grepl(BLACKLIST, rep_df$target, ignore.case=TRUE, perl=TRUE)
  rep_named <- rep[keep]
  rep_named <- rep_named[!duplicated(as.data.frame(rep_named)$target)]
  message("  After blacklist filter: ", length(rep_named), " motifs")

  # Save filtered CSV
  csv_out <- file.path(OUT_DIR,
    paste0("pwmenrich_top_motifs_", cfg$name, "_5k_v2.csv"))
  write.csv(as.data.frame(rep_named), csv_out, row.names=FALSE)
  message("  Saved CSV: ", basename(csv_out))

  # Native plot -- top 20
  N20     <- min(20, length(rep_named))
  pdf_out <- file.path(OUT_DIR,
    paste0("pwmenrich_report_", cfg$name, "_5k_v2.pdf"))
  cairo_pdf(pdf_out, width=12, height=max(10, N20 * 0.65))
  plot(rep_named[1:N20])
  dev.off()
  message("  Saved native plot (top 20): ", basename(pdf_out))

  # Native plot -- top 10
  N10      <- min(10, length(rep_named))
  pdf_top10 <- file.path(OUT_DIR,
    paste0("pwmenrich_top10_", cfg$name, "_5k_v2.pdf"))
  cairo_pdf(pdf_top10, width=11, height=max(8, N10 * 0.70))
  plot(rep_named[1:N10])
  dev.off()
  message("  Saved native plot (top 10): ", basename(pdf_top10))

  # Clean up large objects between contrasts
  rm(gr, seqs, res, rep, rep_named)
  gc()
}

message("\nAll contrasts done. Outputs in: ", OUT_DIR)
