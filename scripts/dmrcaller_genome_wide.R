#!/usr/bin/env Rscript
# =============================================================================
# dmrcaller_genome_wide.R
# Genome-wide DMR calling — SMA Epigenomics Pipeline
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics
# Supervisor: Dr Radu Zabet
#
# Locked parameters (confirmed with Radu, May 2026):
#   Method:     bins
#   Bin size:   300 bp
#   minDiff:    0.20
#   p-value:    0.01
#   min CpGs:   4
#   Permut.:    20 seeds
#   Context:    CpG only
#   Pooling:    poolMethylationDatasets() per condition
#
# Three contrasts:
#   1. ASO_VPA    vs Scramble_CTRL  (combination vs baseline — main result)
#   2. Scramble_VPA vs Scramble_CTRL  (VPA effect alone)
#   3. ASO_CTRL   vs Scramble_CTRL  (ASO effect alone — negative control)
#
# Usage:
#   Rscript scripts/dmrcaller_genome_wide.R
#   (run from project root on Apocrita)
# =============================================================================

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

.libPaths(c("~/R/library", .libPaths()))

# ── Configuration ─────────────────────────────────────────────────────────────
BY_CHR_DIR <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- c(
  "ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3",
  "ASO_VPA_1",  "ASO_VPA_2",  "ASO_VPA_3",
  "Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
  "Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3"
)

CONDITIONS <- list(
  ASO_CTRL     = paste0("ASO_CTRL_",     1:3),
  ASO_VPA      = paste0("ASO_VPA_",      1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)

# Standard chromosomes only
CHROMS <- paste0("chr", c(1:22, "X", "Y"))

# Locked DMR parameters
DMR_PARAMS <- list(
  method       = "bins",
  binSize      = 300,
  minProportionDifference = 0.20,
  pValueThreshold         = 0.01,
  minCytosinesCount       = 4,
  minReadsPerCytosine     = 1,
  context      = "CG",
  permutations = 20
)

# Contrasts: list(name, condition_A, condition_B, description)
CONTRASTS <- list(
  list(
    name    = "ASO_VPA_vs_Scramble_CTRL",
    cond_a  = "ASO_VPA",
    cond_b  = "Scramble_CTRL",
    desc    = "Combination treatment vs baseline — main pleiotropic signal"
  ),
  list(
    name    = "Scramble_VPA_vs_Scramble_CTRL",
    cond_a  = "Scramble_VPA",
    cond_b  = "Scramble_CTRL",
    desc    = "VPA effect alone — expected to overlap with combination"
  ),
  list(
    name    = "ASO_CTRL_vs_Scramble_CTRL",
    cond_a  = "ASO_CTRL",
    cond_b  = "Scramble_CTRL",
    desc    = "ASO alone — negative control, expect minimal DMRs"
  )
)

# ── Step 1: Load per-replicate CpG reports ────────────────────────────────────
message("=== Loading CpG reports ===")
meth_raw <- list()

for (sample in SAMPLES) {
  message("  Loading: ", sample)
  chr_files <- file.path(BY_CHR_DIR,
    paste0(sample, "_", CHROMS, ".CpG_report.txt.gz"))
  chr_files <- chr_files[file.exists(chr_files)]

  if (length(chr_files) == 0) {
    stop("No chr files found for sample: ", sample,
         "\n  Expected in: ", BY_CHR_DIR)
  }

  meth_list <- lapply(chr_files, readBismark)
  meth_raw[[sample]] <- do.call(c, meth_list)
}
message("All samples loaded.")

# ── Step 2: Pool replicates per condition ─────────────────────────────────────
message("\n=== Pooling replicates per condition ===")
meth_pooled <- list()

for (cond in names(CONDITIONS)) {
  message("  Pooling: ", cond)
  reps  <- CONDITIONS[[cond]]
  glist <- GRangesList(lapply(reps, function(s) meth_raw[[s]]))
  meth_pooled[[cond]] <- poolMethylationDatasets(glist)
  message("    CpGs after pooling: ",
          format(length(meth_pooled[[cond]]), big.mark = ","))
}

# ── Step 3: Run DMR calling for each contrast ─────────────────────────────────
message("\n=== Running DMR calling ===")
message("Parameters:")
message("  method:  ", DMR_PARAMS$method)
message("  binSize: ", DMR_PARAMS$binSize, " bp")
message("  minDiff: ", DMR_PARAMS$minProportionDifference)
message("  p-value: ", DMR_PARAMS$pValueThreshold)
message("  minCpGs: ", DMR_PARAMS$minCytosinesCount)
message("  permut.: ", DMR_PARAMS$permutations)
message("  context: ", DMR_PARAMS$context)

results_all <- list()
summary_rows <- list()

for (ct in CONTRASTS) {
  message("\n--- Contrast: ", ct$name, " ---")
  message("    ", ct$desc)

  set.seed(42)
  dmrs <- computeDMRs(
    methylationData1     = meth_pooled[[ct$cond_a]],
    methylationData2     = meth_pooled[[ct$cond_b]],
    regions              = NULL,
    context              = DMR_PARAMS$context,
    method               = DMR_PARAMS$method,
    binSize              = DMR_PARAMS$binSize,
    minProportionDifference = DMR_PARAMS$minProportionDifference,
    pValueThreshold      = DMR_PARAMS$pValueThreshold,
    minCytosinesCount    = DMR_PARAMS$minCytosinesCount,
    minReadsPerCytosine  = DMR_PARAMS$minReadsPerCytosine,
    test                 = "chi-squared",
    pseudocountM         = 1,
    pseudocountN         = 2,
    cores                = 1,
    permutations         = DMR_PARAMS$permutations
  )

  results_all[[ct$name]] <- dmrs
  n_total <- length(dmrs)
  n_hyper <- if (n_total > 0) sum(dmrs$regionType == "gain", na.rm = TRUE) else 0
  n_hypo  <- if (n_total > 0) sum(dmrs$regionType == "loss", na.rm = TRUE) else 0

  message("  DMRs found: ", n_total,
          " (", n_hyper, " hyper, ", n_hypo, " hypo)")

  # Save per-contrast RDS
  rds_path <- file.path(OUT_DIR, paste0("dmr_", ct$name, ".rds"))
  saveRDS(dmrs, rds_path)
  message("  Saved: ", rds_path)

  # Save per-contrast BED
  if (n_total > 0) {
    bed_path <- file.path(OUT_DIR, paste0("dmr_", ct$name, ".bed"))
    bed_df <- data.frame(
      chr        = as.character(seqnames(dmrs)),
      start      = start(dmrs) - 1,
      end        = end(dmrs),
      name       = ct$name,
      score      = round(-log10(dmrs$pValue + 1e-300), 2),
      strand      = ".",
      methDiff   = round(dmrs$regionType == "gain" - dmrs$regionType == "loss", 4),
      type       = dmrs$regionType,
      nCpG       = dmrs$cytosinesCount,
      pValue     = dmrs$pValue
    )
    write.table(bed_df, bed_path, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    message("  Saved: ", bed_path)
  }

  summary_rows[[ct$name]] <- data.frame(
    contrast    = ct$name,
    description = ct$desc,
    n_DMRs      = n_total,
    n_hyper     = n_hyper,
    n_hypo      = n_hypo,
    pct_hypo    = if (n_total > 0) round(100 * n_hypo / n_total, 1) else NA
  )
}

# ── Step 4: Summary table ─────────────────────────────────────────────────────
message("\n=== Summary ===")
summary_df <- do.call(rbind, summary_rows)
print(summary_df)

tsv_path <- file.path(OUT_DIR, "dmr_summary.tsv")
write.table(summary_df, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Summary saved: ", tsv_path)

# ── Step 5: Save all results together ─────────────────────────────────────────
saveRDS(results_all, file.path(OUT_DIR, "dmr_all_contrasts.rds"))
message("\n=== DMR calling complete ===")
message("Outputs in: ", OUT_DIR)
message("Parameters used:")
message("  method=", DMR_PARAMS$method,
        " | binSize=", DMR_PARAMS$binSize,
        " | minDiff=", DMR_PARAMS$minProportionDifference,
        " | p<=", DMR_PARAMS$pValueThreshold,
        " | minCpGs=", DMR_PARAMS$minCytosinesCount,
        " | permutations=", DMR_PARAMS$permutations)
