#!/usr/bin/env Rscript
# ---
# dmrcaller_genome_wide.R
# Genome-wide DMR calling (single-job version) — SMA Epigenomics Pipeline
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics
# Supervisor: Dr Radu Zabet
#
# NOTE: This script was the original genome-wide approach. It is kept for
# reference but the per-chromosome approach (dmrcaller_by_chr.R +
# dmrcaller_combine_chr.R) is now preferred because it avoids the iterative
# merge hang on VPA contrasts and is 72x faster via SLURM parallelisation.
# Use this script only if you need a single-job run for a specific contrast.
#
# Parameters locked with Radu Zabet, 5 May 2026.
#
# Usage:
#   Rscript scripts/dmrcaller_genome_wide.R                     # all 3 contrasts
#   Rscript scripts/dmrcaller_genome_wide.R ASO_VPA_vs_Scramble_CTRL  # single contrast
# ---
 
suppressPackageStartupMessages({
  library(DMRcaller)      # Catoni et al. 2018, NAR
  library(GenomicRanges)  # Bioconductor genomic interval operations
  library(R.utils)        # withTimeout() — prevents infinite merge hangs
})
 
.libPaths(c("~/R/library", .libPaths()))
 
# Paths
# Uses UNMASKED alignment data. Masking only affects chr5:70,924,941-70,953,015
# (SMN1 locus). For all other chromosomes masked/unmasked data are identical.
# Masked data is used in dmrcaller_smn_locus_masked.R for SMN2 locus only.
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
 
# Exclude unplaced contigs, patches and decoys — very low coverage,
# inflate DMR counts with unreliable calls
CHROMS <- paste0("chr", c(1:22, "X", "Y"))
 
# Locked DMR parameters
# All confirmed with Radu Zabet on 5 May 2026.
# Do not change without re-benchmarking and confirming with Radu.
DMR_PARAMS <- list(
 
  # "bins" divides the genome into fixed 300bp windows and tests each for
  # differential methylation. Alternative is "noise_filter" (variable-width).
  # We chose bins because: (1) each window is the same size = interpretable;
  # (2) doesn't require pre-specifying CpG-dense regions;
  # (3) performs well at ~27x pooled coverage.
  method  = "bins",
  binSize = 300,  # benchmarked against chr1 permutation null (benchmark_nb_mingap_radu.R)
 
  # 20% absolute methylation difference threshold. Filters biologically trivial
  # fluctuations. With mean methylation ~80% at SMN, a 20% change is substantial.
  minProportionDifference = 0.20,
 
  pValueThreshold = 0.01,  # standard genome-wide threshold
 
  # Minimum 4 CpGs per 300bp bin. Prevents calling DMRs from single-CpG noise.
  # At average hg38 CpG density (~1 per 107bp), most 300bp bins have 2-3 CpGs;
  # requiring 4 ensures genuine signal.
  minCytosinesCount = 4,
 
  # Minimum 4 reads per CpG position. Confirmed from benchmark_nb_mingap_radu.R.
  # Using the default (1) produced >500,000 spurious DMRs from low-coverage noise.
  # This is the most impactful filter for clean results.
  minReadsPerCytosine = 4,
 
  context = "CG",  # CpG only. CHG/CHH is near-zero in human somatic cells.
 
  # Timeout per contrast in seconds. computeDMRs hangs indefinitely in the
  # iterative merge step when VPA creates thousands of adjacent significant bins.
  # 86400s = 24h. If this fires, use the per-chromosome approach instead.
  timeout_seconds = 86400
)
 
# Contrasts
# Each contrast isolates one experimental variable by holding the other fixed.
# methylationData1 = treatment (cond_a), methylationData2 = reference (cond_b).
CONTRASTS <- list(
  list(
    name   = "ASO_VPA_vs_Scramble_CTRL",
    cond_a = "ASO_VPA",       # nusinersen + VPA (both treatments)
    cond_b = "Scramble_CTRL", # neither treatment — true baseline
    desc   = "Combination vs baseline — primary result"
  ),
  list(
    name   = "Scramble_VPA_vs_Scramble_CTRL",
    cond_a = "Scramble_VPA",  # VPA only, no ASO
    cond_b = "Scramble_CTRL",
    desc   = "VPA alone vs baseline — isolates HDAC inhibitor effect"
  ),
  list(
    name   = "ASO_CTRL_vs_Scramble_CTRL",
    cond_a = "ASO_CTRL",      # nusinersen only, no VPA
    cond_b = "Scramble_CTRL",
    desc   = "Nusinersen alone vs baseline — negative control, expect minimal DMRs"
  )
)
 
# Optional: run single contrast from command line
# Allows each contrast to be submitted as a separate SLURM job so a hang
# in one doesn't block the others. Loads only the 6 samples needed.
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
 
# Step 1: Load per-sample CpG reports
# readBismark() reads bismark --cytosine_report format.
# Critical: this format includes ALL cytosine positions in the genome including
# zero-coverage ones. DMRcaller needs zero-coverage bins for correct statistics.
# If you use a coverage-only format (bedGraph), you will get biased results.
message("=== Loading CpG reports ===")
meth_raw <- list()
 
# Only load conditions needed for selected contrasts — saves ~30min if running
# a single contrast rather than loading all 12 samples unnecessarily
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
  # Load each chromosome separately then concatenate.
  # Loading all chromosomes at once would require loading the full ~28M CpG
  # genome-wide report into RAM in one shot — per-chr keeps memory manageable.
  meth_list          <- lapply(chr_files, readBismark)
  meth_raw[[sample]] <- do.call(c, meth_list)  # concatenate GRanges across chromosomes
}
message("All samples loaded.")
 
# Step 2: Pool replicates
# poolMethylationDatasets() sums methylated + unmethylated read counts across
# all 3 replicates at each CpG position, giving ~27x effective depth.
#
# WHY pool rather than computeDMRsReplicates()?
# computeDMRsReplicates() uses a mixed-effects model which requires sufficient
# per-replicate coverage (ideally >20x) to estimate between-replicate variance.
# Our per-replicate coverage is ~9x — too low for reliable variance estimation.
# Pooling collapses replicates into a single deeper-coverage pseudo-sample;
# confirmed as the appropriate strategy with Radu on 5 May 2026.
message("\n=== Pooling replicates per condition ===")
meth_pooled <- list()
 
for (cond in needed_conds) {
  message("  Pooling: ", cond)
  reps  <- CONDITIONS[[cond]]
  glist <- GRangesList(lapply(reps, function(s) meth_raw[[s]]))
  meth_pooled[[cond]] <- poolMethylationDatasets(glist)
  message("    CpGs after pooling: ",
          format(length(meth_pooled[[cond]]), big.mark = ","))
}
 
# Step 3: DMR calling
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
 
  # Per-contrast checkpoint — if RDS already exists, load and skip recomputation.
  # Essential for rerunning after a timeout: completed contrasts are preserved.
  rds_path <- file.path(OUT_DIR, paste0("dmr_", ct$name, ".rds"))
  if (file.exists(rds_path)) {
    message("  Already exists, loading from: ", rds_path)
    dmrs    <- readRDS(rds_path)
    results_all[[ct$name]] <- dmrs
    n_total <- length(dmrs)
    n_hyper <- if (n_total > 0) sum(dmrs$regionType == "loss", na.rm = TRUE) else 0
    n_hypo  <- if (n_total > 0) sum(dmrs$regionType == "gain", na.rm = TRUE) else 0
    message("  DMRs loaded: ", n_total, " (", n_hyper, " hyper / ", n_hypo, " hypo)")
    summary_rows[[ct$name]] <- data.frame(
      contrast = ct$name, desc = ct$desc, n_DMRs = n_total,
      n_hyper = n_hyper, n_hypo = n_hypo,
      pct_hypo = if (n_total > 0) round(100 * n_hypo / n_total, 1) else NA
    )
    next
  }
 
  # Register MulticoreParam BEFORE calling computeDMRs.
  # DMRcaller uses BiocParallel internally for parallelisation.
  # Default backend is SnowParam (socket-based PSOCK clusters) which
  # causes random hangs on Apocrita SLURM nodes due to port conflicts.
  # MulticoreParam uses fork-based parallelism which is stable on Linux HPC.
  # Must be registered fresh per contrast — not just once at script start.
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))
 
  set.seed(42)  # reproducibility for any stochastic steps inside DMRcaller
 
  # withTimeout() from R.utils wraps computeDMRs with a wall-clock timeout.
  # computeDMRs can hang indefinitely in the iterative merge step when VPA
  # causes thousands of adjacent significant bins. The timeout allows the
  # script to skip a stuck contrast and continue with the next one, rather
  # than blocking for days. If this fires, switch to dmrcaller_by_chr.R.
  dmrs <- tryCatch({
    withTimeout(
      computeDMRs(
        methylationData1        = meth_pooled[[ct$cond_a]],
        methylationData2        = meth_pooled[[ct$cond_b]],
        regions                 = NULL,      # genome-wide, no pre-filtering
        context                 = DMR_PARAMS$context,
        method                  = DMR_PARAMS$method,
        binSize                 = DMR_PARAMS$binSize,
        minProportionDifference = DMR_PARAMS$minProportionDifference,
        pValueThreshold         = DMR_PARAMS$pValueThreshold,
        minCytosinesCount       = DMR_PARAMS$minCytosinesCount,
        minReadsPerCytosine     = DMR_PARAMS$minReadsPerCytosine,
        # score test = Rao test, a large-sample approximation to the LRT.
        # Appropriate at ~27x pooled coverage. Faster than Fisher's exact.
        # "chi-squared" is NOT a valid option in DMRcaller v0.25.1 — will error.
        test  = "score",
        # minGap = 300: minimum gap between DMRs = one full bin width.
        # Prevents the iterative merge from chaining adjacent bins into
        # one enormous DMR spanning megabases. PENDING confirmation from Radu.
        minGap = 300,
        # cores = 4: number of parallel cores for computeDMRs.
        # Uses MulticoreParam registered above.
        cores = 4
      ),
      timeout   = DMR_PARAMS$timeout_seconds,
      onTimeout = "error"
    )
  }, TimeoutException = function(e) {
    message("  TIMEOUT after ", DMR_PARAMS$timeout_seconds, "s — contrast: ", ct$name)
    message("  Switch to per-chromosome approach:")
    message("  bash scripts/submit_dmr_by_chr.sh")
    return(NULL)
  }, error = function(e) {
    message("  ERROR in computeDMRs for ", ct$name, ": ", conditionMessage(e))
    return(NULL)
  })
 
  if (is.null(dmrs)) next  # timeout or error — move to next contrast
 
  results_all[[ct$name]] <- dmrs
  n_total <- length(dmrs)
 
  # regionType convention: "gain" = hypomethylated in treatment (proportion1 < proportion2)
  #                        "loss" = hypermethylated in treatment (proportion1 > proportion2)
  # Confirmed by inspecting chr13 DMRs: proportion1 (ASO_VPA) < proportion2 (Scramble_CTRL)
  # yet regionType = "gain" — so "gain" means the treatment GAINED hypomethylation vs reference
  n_hyper <- if (n_total > 0) sum(dmrs$regionType == "loss", na.rm = TRUE) else 0
  n_hypo  <- if (n_total > 0) sum(dmrs$regionType == "gain", na.rm = TRUE) else 0
  message("  DMRs: ", n_total, " (", n_hyper, " hyper / ", n_hypo, " hypo)")
 
  # Save RDS immediately after each contrast completes — don't wait for all 3.
  # If the next contrast hangs, this one is already persisted to disk.
  saveRDS(dmrs, rds_path)
  message("  RDS saved: ", rds_path)
 
  # BED file for IGV / genome browser visualisation
  # start - 1: converts from 1-based GRanges to 0-based BED coordinates
  if (n_total > 0) {
    bed_path <- file.path(OUT_DIR, paste0("dmr_", ct$name, ".bed"))
    bed_df <- data.frame(
      chr      = as.character(seqnames(dmrs)),
      start    = start(dmrs) - 1,
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
    message("  BED saved: ", bed_path)
  }
 
  summary_rows[[ct$name]] <- data.frame(
    contrast = ct$name, desc = ct$desc, n_DMRs = n_total,
    n_hyper  = n_hyper, n_hypo = n_hypo,
    pct_hypo = if (n_total > 0) round(100 * n_hypo / n_total, 1) else NA
  )
}
 
# Summary
if (length(summary_rows) > 0) {
  message("\n=== Summary ===")
  summary_df <- do.call(rbind, summary_rows)
  print(summary_df)
  write.table(summary_df, file.path(OUT_DIR, "dmr_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Summary saved: ", file.path(OUT_DIR, "dmr_summary.tsv"))
}
 
# Save combined RDS with all three contrasts — convenient for quick
# interactive exploration in R without loading 3 separate files
if (length(results_all) > 0) {
  saveRDS(results_all, file.path(OUT_DIR, "dmr_all_contrasts.rds"))
}
 
message("\nDone. Outputs in: ", OUT_DIR)
