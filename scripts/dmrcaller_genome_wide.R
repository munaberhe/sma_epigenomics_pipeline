#!/usr/bin/env Rscript
# =============================================================================
# dmrcaller_genome_wide.R
# Genome-wide DMR calling (single-job version) — SMA Epigenomics Pipeline
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics
# Supervisor: Dr Radu Zabet
#
# NOTE: kept for reference. The per-chromosome approach (dmrcaller_by_chr.R +
# dmrcaller_combine_chr.R) is now preferred — it avoids the iterative merge
# hang on VPA contrasts and is 72x faster via SLURM parallelisation.
# Use this script for a quick single-contrast run or for debugging.
#
# Parameters locked with Radu Zabet, 5 May 2026.
#
# Usage:
#   Rscript scripts/dmrcaller_genome_wide.R                          # all 3
#   Rscript scripts/dmrcaller_genome_wide.R ASO_VPA_vs_Scramble_CTRL # single
# =============================================================================
 
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(R.utils)   # withTimeout() — kills computeDMRs if merge hangs
})
 
.libPaths(c("~/R/library", .libPaths()))
 
# Unmasked alignment data. Masking only affects chr5:70,924,941-70,953,015.
# For all other chromosomes masked and unmasked data are identical.
# Masked data is used separately in dmrcaller_smn_locus_masked.R.
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
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = paste0("ASO_VPA_",       1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)
 
# Exclude unplaced contigs and patches — low coverage, inflate DMR counts
CHROMS <- paste0("chr", c(1:22, "X", "Y"))
 
# ── Locked DMR parameters (confirmed with Radu, 5 May 2026) ──────────────────
DMR_PARAMS <- list(
  # "bins": fixed 300bp windows across the genome. Chosen over "noise_filter"
  # because every window is the same size (interpretable), doesn't require
  # pre-specifying CpG-dense regions, and performs well at ~27x pooled depth.
  method  = "bins",
  binSize = 300,   # benchmarked on chr1 permutation null (benchmark_nb_mingap_radu.R)
 
  minProportionDifference = 0.20,   # 20% absolute difference; filters trivial changes
  pValueThreshold         = 0.01,
 
  # min 4 CpGs per 300bp bin. At hg38 average density (~1 CpG per 107bp),
  # most 300bp bins have 2-3 CpGs; requiring 4 filters the noisiest bins.
  minCytosinesCount = 4,
 
  # min 4 reads per CpG. Default of 1 produced >500k spurious DMRs from
  # low-coverage noise. Confirmed from benchmark_nb_mingap_radu.R.
  minReadsPerCytosine = 4,
 
  context = "CG",   # CpG only; CHG/CHH near-zero in human somatic cells
 
  # 86400s = 24h timeout. If computeDMRs hasn't returned by then it's stuck
  # in the iterative merge. Switch to dmrcaller_by_chr.R instead.
  timeout_seconds = 86400
)
 
# ── Contrasts ─────────────────────────────────────────────────────────────────
CONTRASTS <- list(
  list(
    name   = "ASO_VPA_vs_Scramble_CTRL",
    cond_a = "ASO_VPA",
    cond_b = "Scramble_CTRL",
    desc   = "Combination vs baseline — primary result"
  ),
  list(
    name   = "Scramble_VPA_vs_Scramble_CTRL",
    cond_a = "Scramble_VPA",
    cond_b = "Scramble_CTRL",
    desc   = "VPA alone vs baseline — isolates HDAC inhibitor effect"
  ),
  list(
    name   = "ASO_CTRL_vs_Scramble_CTRL",
    cond_a = "ASO_CTRL",
    cond_b = "Scramble_CTRL",
    desc   = "Nusinersen alone vs baseline — negative control"
  )
)
 
# Single contrast mode — loads only the 6 samples needed, saves ~30min
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  contrast_filter <- args[1]
  CONTRASTS <- Filter(function(ct) ct$name == contrast_filter, CONTRASTS)
  if (length(CONTRASTS) == 0) {
    stop("Unknown contrast: ", contrast_filter,
         "\nValid: ASO_VPA_vs_Scramble_CTRL, ",
         "Scramble_VPA_vs_Scramble_CTRL, ASO_CTRL_vs_Scramble_CTRL")
  }
  message("Running single contrast: ", contrast_filter)
}
 
# ── Step 1: Load CpG reports ──────────────────────────────────────────────────
# readBismark() reads bismark --cytosine_report format. Includes zero-coverage
# positions — DMRcaller needs these for correct bin statistics.
message("=== Loading CpG reports ===")
meth_raw <- list()
 
needed_conds   <- unique(unlist(lapply(CONTRASTS, function(ct) c(ct$cond_a, ct$cond_b))))
needed_samples <- unique(unlist(CONDITIONS[needed_conds]))
 
for (sample in needed_samples) {
  message("  Loading: ", sample)
  chr_files <- file.path(BY_CHR_DIR,
    paste0(sample, "_", CHROMS, ".CpG_report.txt.gz"))
  chr_files <- chr_files[file.exists(chr_files)]
  if (length(chr_files) == 0) {
    stop("No chr files found for sample: ", sample,
         "\n  Run 04_split_by_chr.sh first.")
  }
  meth_list          <- lapply(chr_files, readBismark)
  meth_raw[[sample]] <- do.call(c, meth_list)
}
message("All samples loaded.")
 
# ── Step 2: Pool replicates ───────────────────────────────────────────────────
# poolMethylationDatasets() sums reads across 3 replicates → ~27x depth.
# computeDMRsReplicates() not used: it needs >20x per replicate for reliable
# variance estimation; we only have ~9x.
message("\n=== Pooling replicates per condition ===")
meth_pooled <- list()
 
for (cond in needed_conds) {
  message("  Pooling: ", cond)
  reps  <- CONDITIONS[[cond]]
  glist <- GRangesList(lapply(reps, function(s) meth_raw[[s]]))
  meth_pooled[[cond]] <- poolMethylationDatasets(glist)
  message("    CpGs: ", format(length(meth_pooled[[cond]]), big.mark = ","))
}
 
# ── Step 3: DMR calling ───────────────────────────────────────────────────────
message("\n=== Running DMR calling ===")
message("  method=", DMR_PARAMS$method,
        " | binSize=", DMR_PARAMS$binSize,
        " | minDiff=", DMR_PARAMS$minProportionDifference,
        " | p<=", DMR_PARAMS$pValueThreshold,
        " | minCpGs=", DMR_PARAMS$minCytosinesCount,
        " | minReads=", DMR_PARAMS$minReadsPerCytosine,
        " | timeout=", DMR_PARAMS$timeout_seconds, "s")
 
results_all  <- list()
summary_rows <- list()
 
for (ct in CONTRASTS) {
  message("\n--- Contrast: ", ct$name, " ---")
  message("    ", ct$desc)
 
  # Checkpoint — skip if RDS already exists. Completed contrasts survive reruns.
  rds_path <- file.path(OUT_DIR, paste0("dmr_", ct$name, ".rds"))
  if (file.exists(rds_path)) {
    message("  Already exists, loading: ", rds_path)
    dmrs    <- readRDS(rds_path)
    results_all[[ct$name]] <- dmrs
    n_total <- length(dmrs)
    n_hyper <- if (n_total > 0) sum(dmrs$regionType == "loss", na.rm = TRUE) else 0
    n_hypo  <- if (n_total > 0) sum(dmrs$regionType == "gain", na.rm = TRUE) else 0
    message("  DMRs: ", n_total, " (", n_hyper, " hyper / ", n_hypo, " hypo)")
    summary_rows[[ct$name]] <- data.frame(
      contrast = ct$name, desc = ct$desc, n_DMRs = n_total,
      n_hyper = n_hyper, n_hypo = n_hypo,
      pct_hypo = if (n_total > 0) round(100 * n_hypo / n_total, 1) else NA
    )
    next
  }
 
  # MulticoreParam (fork-based) instead of default SnowParam (socket-based).
  # SnowParam causes random hangs on Apocrita SLURM nodes due to port conflicts.
  # Must be registered per-contrast, not once at script start.
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))
 
  set.seed(42)
 
  # withTimeout() wraps computeDMRs with a wall-clock limit.
  # VPA causes so many adjacent significant bins that the iterative merge
  # never converges on a genome-wide run. The timeout lets the script skip
  # a stuck contrast rather than blocking for days.
  dmrs <- tryCatch({
    withTimeout(
      computeDMRs(
        methylationData1        = meth_pooled[[ct$cond_a]],
        methylationData2        = meth_pooled[[ct$cond_b]],
        regions                 = NULL,
        context                 = DMR_PARAMS$context,
        method                  = DMR_PARAMS$method,
        binSize                 = DMR_PARAMS$binSize,
        minProportionDifference = DMR_PARAMS$minProportionDifference,
        pValueThreshold         = DMR_PARAMS$pValueThreshold,
        minCytosinesCount       = DMR_PARAMS$minCytosinesCount,
        minReadsPerCytosine     = DMR_PARAMS$minReadsPerCytosine,
        # score test = Rao test; "chi-squared" is not valid in v0.25.1
        test   = "score",
        # minGap = 1 bin width; prevents chaining adjacent bins into
        # megabase-spanning DMRs on VPA contrasts — PENDING Radu confirmation
        minGap = 300,
        cores  = 4
      ),
      timeout   = DMR_PARAMS$timeout_seconds,
      onTimeout = "error"
    )
  }, TimeoutException = function(e) {
    message("  TIMEOUT — switch to per-chromosome approach:")
    message("  bash scripts/submit_dmr_by_chr.sh")
    return(NULL)
  }, error = function(e) {
    message("  ERROR: ", conditionMessage(e))
    return(NULL)
  })
 
  if (is.null(dmrs)) next
 
  results_all[[ct$name]] <- dmrs
  n_total <- length(dmrs)
  # "gain" = hypomethylated in treatment, "loss" = hypermethylated
  n_hyper <- if (n_total > 0) sum(dmrs$regionType == "loss", na.rm = TRUE) else 0
  n_hypo  <- if (n_total > 0) sum(dmrs$regionType == "gain", na.rm = TRUE) else 0
  message("  DMRs: ", n_total, " (", n_hyper, " hyper / ", n_hypo, " hypo)")
 
  # Save per-contrast immediately — don't wait for all 3 to finish
  saveRDS(dmrs, rds_path)
  message("  RDS: ", rds_path)
 
  if (n_total > 0) {
    bed_path <- file.path(OUT_DIR, paste0("dmr_", ct$name, ".bed"))
    bed_df <- data.frame(
      chr      = as.character(seqnames(dmrs)),
      start    = start(dmrs) - 1,   # 1-based GRanges → 0-based BED
      end      = end(dmrs),
      name     = ct$name,
      score    = round(-log10(dmrs$pValue + 1e-300), 2),
      strand   = ".",
      methDiff = round(as.integer(dmrs$regionType == "loss") -
                       as.integer(dmrs$regionType == "gain"), 4),
      type     = dmrs$regionType,
      nCpG     = dmrs$cytosinesCount,
      pValue   = dmrs$pValue
    )
    write.table(bed_df, bed_path, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    message("  BED: ", bed_path)
  }
 
  summary_rows[[ct$name]] <- data.frame(
    contrast = ct$name, desc = ct$desc, n_DMRs = n_total,
    n_hyper = n_hyper, n_hypo = n_hypo,
    pct_hypo = if (n_total > 0) round(100 * n_hypo / n_total, 1) else NA
  )
}
 
if (length(summary_rows) > 0) {
  message("\n=== Summary ===")
  summary_df <- do.call(rbind, summary_rows)
  print(summary_df)
  write.table(summary_df, file.path(OUT_DIR, "dmr_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
 
if (length(results_all) > 0) {
  saveRDS(results_all, file.path(OUT_DIR, "dmr_all_contrasts.rds"))
}
 
message("\nDone. Outputs in: ", OUT_DIR)
