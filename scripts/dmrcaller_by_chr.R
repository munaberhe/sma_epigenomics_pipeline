#!/usr/bin/env Rscript
# =============================================================================
# dmrcaller_by_chr.R
# Per-chromosome DMR calling — SMA Epigenomics Pipeline
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics
# Supervisor: Dr Radu Zabet
#
# WHY per-chromosome instead of genome-wide in one shot?
# computeDMRs() runs an iterative merge step after calling DMRs on each
# chromosome. When VPA causes widespread hypomethylation (hundreds of
# thousands of adjacent significant bins), this merge hangs indefinitely.
# Splitting per-chromosome gives computeDMRs() a much smaller problem —
# finishes in minutes per chromosome instead of hanging for days.
# Results are combined by dmrcaller_combine_chr.R.
#
# Usage: Rscript scripts/dmrcaller_by_chr.R <contrast> <chr>
# Example: Rscript scripts/dmrcaller_by_chr.R ASO_VPA_vs_Scramble_CTRL chr13
#
# Submitted as 72 parallel SLURM jobs (24 chr x 3 contrasts) via
# scripts/submit_dmr_by_chr.sh
# =============================================================================
 
suppressPackageStartupMessages({
  library(DMRcaller)      # Catoni et al. 2018, NAR
  library(GenomicRanges)  # Bioconductor genomic interval operations
})
 
.libPaths(c("~/R/library", .libPaths()))
 
# ── Parse args ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript dmrcaller_by_chr.R <contrast> <chr>")
 
CONTRAST <- args[1]  # e.g. "ASO_VPA_vs_Scramble_CTRL"
CHR      <- args[2]  # e.g. "chr13"
 
BY_CHR_DIR <- "results/alignments/bs/by_chr"
OUT_DIR    <- file.path("results/dmr/by_chr", CONTRAST)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
 
# ── Checkpoint ────────────────────────────────────────────────────────────────
# Skip if already done — means you can resubmit the full 72-job array after
# a failure and only the missing chromosomes will actually run.
rds_path <- file.path(OUT_DIR, paste0("dmr_", CONTRAST, "_", CHR, ".rds"))
if (file.exists(rds_path)) {
  message("Already done: ", rds_path)
  quit(status = 0)
}
 
# ── Sample and contrast definitions ──────────────────────────────────────────
CONDITIONS <- list(
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = paste0("ASO_VPA_",       1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)
 
# methylationData1 = treatment (cond_a), methylationData2 = reference (cond_b)
# regionType "gain" = proportion1 < proportion2 = hypomethylated in treatment
# regionType "loss" = proportion1 > proportion2 = hypermethylated in treatment
# Counterintuitive naming confirmed by checking chr13 proportions manually.
CONTRASTS <- list(
  ASO_VPA_vs_Scramble_CTRL       = c("ASO_VPA",      "Scramble_CTRL"),
  Scramble_VPA_vs_Scramble_CTRL  = c("Scramble_VPA", "Scramble_CTRL"),
  ASO_CTRL_vs_Scramble_CTRL      = c("ASO_CTRL",     "Scramble_CTRL")
)
 
if (!CONTRAST %in% names(CONTRASTS)) stop("Unknown contrast: ", CONTRAST)
cond_a <- CONTRASTS[[CONTRAST]][1]
cond_b <- CONTRASTS[[CONTRAST]][2]
 
# ── Load and pool one chromosome ──────────────────────────────────────────────
# readBismark() reads bismark --cytosine_report format. The key requirement:
# this format includes zero-coverage positions. DMRcaller needs these for
# correct bin statistics — a bedGraph (coverage-only) would give biased results.
#
# poolMethylationDatasets() sums methylated + unmethylated reads across the 3
# replicates at each CpG, giving ~27x depth from ~9x per replicate.
# Using computeDMRsReplicates() instead would require >20x per replicate for
# reliable variance estimation — we don't have that.
load_chr <- function(condition) {
  samples <- CONDITIONS[[condition]]
  grs <- lapply(samples, function(s) {
    f <- file.path(BY_CHR_DIR, paste0(s, "_", CHR, ".CpG_report.txt.gz"))
    if (!file.exists(f)) stop("Missing: ", f, "\nRun 04_split_by_chr.sh first.")
    readBismark(f)
  })
  poolMethylationDatasets(GRangesList(grs))
}
 
message("Loading ", CHR, " for ", CONTRAST)
meth_a <- load_chr(cond_a)
meth_b <- load_chr(cond_b)
message("Pooled CpGs: ", format(length(meth_a), big.mark = ","))
 
# ── DMR calling ───────────────────────────────────────────────────────────────
# computeDMRs() bins method:
#   1. Divide chromosome into 300bp windows
#   2. Score test (Rao) for differential methylation per bin
#   3. Filter by minProportionDifference, pValueThreshold, minCpGs, minReads
#   4. Merge adjacent significant bins (kept tractable here because we're
#      only merging within one chromosome, not across the whole genome)
#
# regions = NULL: test the full chromosome, no pre-filtering
# cores = 1: we're already parallelising at chromosome level (72 SLURM jobs),
#            so single-threaded per job is correct. BiocParallel's SnowParam
#            also causes random hangs on Apocrita SLURM nodes.
set.seed(42)
dmrs <- computeDMRs(
  methylationData1        = meth_a,
  methylationData2        = meth_b,
  regions                 = NULL,
  context                 = "CG",     # CpG only; CHG/CHH near-zero in human somatic cells
  method                  = "bins",   # confirmed with Radu 5 May 2026
  binSize                 = 300,      # benchmarked on chr1 permutation null
  minProportionDifference = 0.20,     # minimum 20% absolute methylation difference
  pValueThreshold         = 0.01,
  minCytosinesCount       = 4,        # minimum CpGs per bin; prevents single-CpG noise
  minReadsPerCytosine     = 4,        # confirmed from benchmark_nb_mingap_radu.R
  minGap                  = 300,      # = 1 bin width; prevents chaining adjacent bins
                                      # into megabase-spanning DMRs — PENDING Radu confirmation
  test                    = "score",  # Rao test; "chi-squared" is not valid in v0.25.1
  cores                   = 1
)
 
message("DMRs on ", CHR, ": ", length(dmrs))
saveRDS(dmrs, rds_path)
message("Saved: ", rds_path)
