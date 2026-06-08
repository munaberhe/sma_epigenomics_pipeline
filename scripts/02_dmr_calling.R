#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

# Per-chromosome DMR calling.
# We run per-chromosome rather than genome-wide because the iterative merge
# step in computeDMRs hangs for 18+ hours on VPA contrasts genome-wide.
# Results are combined by 06b_dmrcaller_combine_chr.R afterwards.
#
# Usage: Rscript 06_dmrcaller_by_chr.R <contrast> <chr>
# Submitted as 72 parallel SLURM jobs via submit_dmr_by_chr.sh

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 06_dmrcaller_by_chr.R <contrast> <chr>")
CONTRAST <- args[1]
CHR      <- args[2]

BY_CHR_DIR <- "results/alignments/bs/by_chr"
OUT_DIR    <- file.path("results/dmr/by_chr", CONTRAST)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# skip if already done — important for partial reruns
rds_path <- file.path(OUT_DIR, paste0("dmr_", CONTRAST, "_", CHR, ".rds"))
if (file.exists(rds_path)) {
  message("Already done: ", rds_path)
  quit(status = 0)
}

CONDITIONS <- list(
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = paste0("ASO_VPA_",       1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)

# contrasts: cond_a = treatment, cond_b = reference
# note: in DMRcaller output, "gain" means cond_a is LOWER than cond_b
# (i.e. hypomethylated in treatment) — counterintuitive naming
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

# load one chromosome for a condition — reads all 3 replicates and pools them
# pooling sums readsM and readsN across reps giving ~27x effective depth
# we pool rather than using computeDMRsReplicates() because ~9x per replicate
# is too low for the mixed-effects model that function requires
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

# DMR calling with locked parameters (benchmarked on chr1 permutation null)
# cores=1 because BiocParallel parallel backend causes random hangs on Apocrita
# SLURM nodes — we parallelise at the chromosome level instead
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
