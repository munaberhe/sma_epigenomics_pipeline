#!/usr/bin/env Rscript
# ---
# dmrcaller_by_chr.R
# Per-chromosome DMR calling — SMA Epigenomics Pipeline
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics
# Supervisor: Dr Radu Zabet
#
# WHY per-chromosome instead of genome-wide in one shot?
# computeDMRs() runs an iterative merge step after calling DMRs on each
# chromosome. When VPA causes widespread hypomethylation (hundreds of
# thousands of adjacent significant bins), this merge hangs indefinitely.
# By splitting per-chromosome we give computeDMRs() a much smaller problem
# to merge — it finishes in minutes per chromosome instead of hanging for days.
# Results are combined afterwards by dmrcaller_combine_chr.R.
#
# Usage: Rscript scripts/dmrcaller_by_chr.R <contrast> <chr>
# Example: Rscript scripts/dmrcaller_by_chr.R ASO_VPA_vs_Scramble_CTRL chr13
#
# Submitted as 72 parallel SLURM jobs (24 chr x 3 contrasts) via
# scripts/submit_dmr_by_chr.sh
# ---
 
suppressPackageStartupMessages({
  library(DMRcaller)      # Catoni et al. 2018, NAR — Radu's lab tool
  library(GenomicRanges)  # Bioconductor genomic interval operations
})
 
# Put our R library first so project packages take precedence over system ones
.libPaths(c("~/R/library", .libPaths()))
 
# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript dmrcaller_by_chr.R <contrast> <chr>")
 
CONTRAST <- args[1]  # e.g. "ASO_VPA_vs_Scramble_CTRL"
CHR      <- args[2]  # e.g. "chr13"
 
# Input: per-chromosome CpG reports from 04_split_by_chr.sh
# These are CpG-only (CHG/CHH filtered out during splitting)
BY_CHR_DIR <- "results/alignments/bs/by_chr"
 
OUT_DIR <- file.path("results/dmr/by_chr", CONTRAST)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
 
# Checkpoint
# If this chromosome's RDS already exists, skip it entirely.
# This is critical for reruns — if 71/72 jobs completed and one failed,
# resubmitting the full array would redo all 71. The checkpoint prevents that.
rds_path <- file.path(OUT_DIR, paste0("dmr_", CONTRAST, "_", CHR, ".rds"))
if (file.exists(rds_path)) {
  message("Already done: ", rds_path)
  quit(status = 0)
}
 
# Sample and contrast definitions
# Each condition has 3 replicates. These are pooled before DMR calling
# (see load_chr() below).
CONDITIONS <- list(
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = paste0("ASO_VPA_",       1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)
 
# Each contrast is a pairwise comparison isolating one experimental variable.
# cond_a = treatment (methylationData1 in computeDMRs)
# cond_b = reference (methylationData2 in computeDMRs)
# IMPORTANT: "gain" in regionType means cond_a is LOWER than cond_b
#            i.e. "gain" = hypomethylated in treatment (counterintuitive naming)
#            "loss" = hypermethylated in treatment
#            This is a DMRcaller convention — confirmed by checking proportion1 vs proportion2
CONTRASTS <- list(
  ASO_VPA_vs_Scramble_CTRL       = c("ASO_VPA",      "Scramble_CTRL"),
  Scramble_VPA_vs_Scramble_CTRL  = c("Scramble_VPA", "Scramble_CTRL"),
  ASO_CTRL_vs_Scramble_CTRL      = c("ASO_CTRL",     "Scramble_CTRL")
)
 
if (!CONTRAST %in% names(CONTRASTS)) stop("Unknown contrast: ", CONTRAST)
cond_a <- CONTRASTS[[CONTRAST]][1]
cond_b <- CONTRASTS[[CONTRAST]][2]
 
# Load and pool one chromosome
# readBismark() reads the genome-wide cytosine report format produced by
# bismark --cytosine_report. Critically, this format includes zero-coverage
# positions — DMRcaller needs these for correct bin statistics. A coverage-only
# file (bedGraph) would give wrong results.
#
# poolMethylationDatasets() sums methylated + unmethylated read counts across
# the 3 replicates at each CpG position. This gives ~27x effective pooled depth
# from ~9x per replicate, which is required for reliable bin-based DMR calling.
# We use pooling rather than computeDMRsReplicates() because per-replicate
# coverage (~9x) is too low for the mixed-effects model that function uses.
load_chr <- function(condition) {
  samples <- CONDITIONS[[condition]]
  grs <- lapply(samples, function(s) {
    f <- file.path(BY_CHR_DIR, paste0(s, "_", CHR, ".CpG_report.txt.gz"))
    if (!file.exists(f)) stop("Missing: ", f,
                               "\nRun 04_split_by_chr.sh first.")
    readBismark(f)
  })
  # GRangesList wraps the 3 individual GRanges objects so pooling can
  # iterate over them — required input format for poolMethylationDatasets()
  poolMethylationDatasets(GRangesList(grs))
}
 
message("Loading ", CHR, " for ", CONTRAST)
meth_a <- load_chr(cond_a)  # treatment condition
meth_b <- load_chr(cond_b)  # reference condition
message("Pooled CpGs: ", format(length(meth_a), big.mark = ","))
 
# DMR calling
# computeDMRs() implements the bins method (Catoni et al. 2018):
#   1. Divide chromosome into fixed 300bp bins
#   2. For each bin: count methylated/total reads in each condition
#   3. Score test (Rao test) for differential methylation
#   4. Filter by minProportionDifference, pValueThreshold, minCytosinesCount,
#      minReadsPerCytosine
#   5. Merge adjacent significant bins (this is the step that hangs genome-wide
#      for VPA — per-chromosome keeps it tractable)
#
# regions = NULL: test the full chromosome, no pre-filtering by region
# cores = 1: BiocParallel's parallel backend causes random hangs on Apocrita
#            SLURM nodes. Since we're already parallelising at the chromosome
#            level (72 SLURM jobs), single-threaded per job is the right choice.
set.seed(42)  # reproducibility for any stochastic steps inside DMRcaller
dmrs <- computeDMRs(
  methylationData1        = meth_a,
  methylationData2        = meth_b,
  regions                 = NULL,
  context                 = "CG",          # CpG only; CHG/CHH near-zero in human somatic cells
  method                  = "bins",        # fixed-width windows; confirmed with Radu 5 May 2026
  binSize                 = 300,           # benchmarked on chr1 permutation null
  minProportionDifference = 0.20,          # minimum 20% absolute methylation difference
  pValueThreshold         = 0.01,          # score test p-value threshold
  minCytosinesCount       = 4,             # minimum CpGs per bin; prevents single-CpG noise
  minReadsPerCytosine     = 4,             # minimum reads per CpG; confirmed from benchmark scripts
  minGap                  = 300,           # minimum gap between DMRs = 1 bin width
                                           # prevents infinite iterative merge on VPA contrasts
                                           # PENDING confirmation from Radu
  test                    = "score",       # Rao score test; appropriate at ~27x pooled depth
                                           # "chi-squared" is not a valid option in this version
  cores                   = 1             # see note above about BiocParallel on Apocrita
)
 
message("DMRs on ", CHR, ": ", length(dmrs))
 
# Save as RDS — primary output format. Read by dmrcaller_combine_chr.R
# which concatenates all per-chromosome results into genome-wide GRanges objects.
saveRDS(dmrs, rds_path)
message("Saved: ", rds_path)
