#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
# Per-chromosome DMR calling.
# FIX (2026-06-19): now also saves the full set of tested windows (before
# p-value filtering) so downstream annotation enrichment (07_grant_fig1cd_*)
# can use the correct background instead of falling back to a genome-wide
# CpG background that doesn't match the actual DMRcaller search space.
#
# Usage: Rscript 02_dmr_calling.R <contrast> <chr>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 02_dmr_calling.R <contrast> <chr>")
CONTRAST <- args[1]
CHR      <- args[2]

BY_CHR_DIR <- "results/alignments/bs/by_chr"
OUT_DIR    <- file.path("results/dmr/by_chr", CONTRAST)
TW_DIR     <- "results/dmr/tested_windows"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TW_DIR,  recursive = TRUE, showWarnings = FALSE)

rds_path <- file.path(OUT_DIR, paste0("dmr_", CONTRAST, "_", CHR, ".rds"))
tw_path  <- file.path(TW_DIR, paste0("tested_windows_", CONTRAST, "_", CHR, ".rds"))

if (file.exists(rds_path) && file.exists(tw_path)) {
  message("Already done: ", rds_path, " and ", tw_path)
  quit(status = 0)
}

CONDITIONS <- list(
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = paste0("ASO_VPA_",       1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)

CONTRASTS <- list(
  ASO_VPA_vs_Scramble_CTRL      = c("ASO_VPA",      "Scramble_CTRL"),
  Scramble_VPA_vs_Scramble_CTRL = c("Scramble_VPA", "Scramble_CTRL"),
  ASO_CTRL_vs_Scramble_CTRL     = c("ASO_CTRL",     "Scramble_CTRL"),
  ASO_VPA_vs_ASO_CTRL           = c("ASO_VPA",      "ASO_CTRL"),
  ASO_VPA_vs_Scramble_VPA       = c("ASO_VPA",      "Scramble_VPA")
)
if (!CONTRAST %in% names(CONTRASTS)) stop("Unknown contrast: ", CONTRAST)
cond_a <- CONTRASTS[[CONTRAST]][1]
cond_b <- CONTRASTS[[CONTRAST]][2]

load_chr <- function(condition) {
  grs <- lapply(CONDITIONS[[condition]], function(s) {
    f <- file.path(BY_CHR_DIR, paste0(s, "_", CHR, ".CpG_report.txt.gz"))
    if (!file.exists(f)) stop("Missing: ", f)
    readBismark(f)
  })
  poolMethylationDatasets(GRangesList(grs))
}

message("Loading ", CHR, " for ", CONTRAST)
meth_a <- load_chr(cond_a)
meth_b <- load_chr(cond_b)
message("Pooled CpGs: ", format(length(meth_a), big.mark=","))

# Build the tested windows independently of computeDMRs' internal filtering.
# These are the 300bp bins that meet the minimum coverage/cytosine criteria
# and were therefore eligible to be tested -- this is the correct background
# for obs/exp annotation enrichment, NOT a genome-wide CpG background.
message("Building tested-window background...")
chr_len <- max(end(meth_a), end(meth_b))
bins <- GRanges(seqnames = CHR,
                ranges = IRanges(start = seq(1, chr_len, by = 300),
                                 width = 300))
# Count cytosines and reads per bin for both conditions, keep bins meeting
# the same minCytosinesCount / minReadsPerCytosine thresholds used in
# computeDMRs below, so the background matches the actual tested space.
ov_a <- findOverlaps(meth_a, bins)
ov_b <- findOverlaps(meth_b, bins)
cyto_per_bin_a <- tabulate(subjectHits(ov_a), nbins = length(bins))
cyto_per_bin_b <- tabulate(subjectHits(ov_b), nbins = length(bins))
reads_a <- sapply(split(meth_a$readsN[queryHits(ov_a)], subjectHits(ov_a)), sum)
reads_b <- sapply(split(meth_b$readsN[queryHits(ov_b)], subjectHits(ov_b)), sum)
mean_reads_a <- rep(0, length(bins)); mean_reads_a[as.integer(names(reads_a))] <-
  reads_a / pmax(cyto_per_bin_a[as.integer(names(reads_a))], 1)
mean_reads_b <- rep(0, length(bins)); mean_reads_b[as.integer(names(reads_b))] <-
  reads_b / pmax(cyto_per_bin_b[as.integer(names(reads_b))], 1)

tested <- bins[cyto_per_bin_a >= 4 & cyto_per_bin_b >= 4 &
               mean_reads_a >= 4 & mean_reads_b >= 4]
message("Tested windows on ", CHR, ": ", length(tested))
saveRDS(tested, tw_path)
message("Saved tested windows: ", tw_path)

# DMR calling parameters
set.seed(42)
dmrs <- computeDMRs(
  methylationData1        = meth_a,
  methylationData2        = meth_b,
  regions                 = NULL,
  context                 = "CG",
  method                  = "bins",
  binSize                 = 300,
  minProportionDifference = 0.20,
  pValueThreshold         = 0.01,
  minCytosinesCount       = 4,
  minReadsPerCytosine     = 4,
  minGap                  = 300,
  test                    = "score",
  cores                   = 1
)
message("DMRs on ", CHR, ": ", length(dmrs))
saveRDS(dmrs, rds_path)
message("Saved: ", rds_path)
