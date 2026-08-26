# labelswap_real_helpers.R
# Real replicate-label permutation for the ASO_VPA vs ASO_CTRL contrast on
# chr1. Replaces the previously broken A<->B swap (which is symmetric under
# the score / Fisher tests used by DMRcaller and so produces n_scrambled
# == n_real exactly, with delta = 0 and ratio = 1 for every cell).
#
# Design
# 6 ASO replicates: ASO_VPA_{1,2,3} + ASO_CTRL_{1,2,3}.
# A 3-vs-3 partition of these 6 replicates gives one (group A pool, group B
# pool) pair. choose(6,3) = 20 partitions exist. Each partition appears as
# its own mirror (A vs B == B vs A under the symmetric DMR test), so there
# are 10 UNIQUE pairings: the identity (real) + 9 non-identity permutations.
#
# We enumerate all 10 exhaustively. No sampling, no bootstrap, no seed
# dependence -- this is the complete null over the experimental design.
#
# Pooling math: readBismarkPool() from DMRcaller, identical to the original
# build_chr1_cache.R. This means the identity permutation (real call)
# reproduces the existing n_real numbers exactly.

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(BiocParallel)
  library(parallel)
})

# paths and constants (match parameter_benchmark_helpers.R)
COV_DIR     <- "results/alignments/bs/by_chr"
OUT_DIR     <- "results/dmr_benchmark_labelswap_real"
CHROM       <- "chr1"
REGION_END  <- 248956422
PERM_DIR    <- file.path(OUT_DIR, "perm_rds")
dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PERM_DIR, recursive = TRUE, showWarnings = FALSE)

ASO_VPA_REPS  <- paste0("ASO_VPA_",  1:3)
ASO_CTRL_REPS <- paste0("ASO_CTRL_", 1:3)
ALL_6         <- c(ASO_VPA_REPS, ASO_CTRL_REPS)  # 1:6 maps to these names

# helpers

rep_path <- function(sample_name) {
  file.path(COV_DIR, paste0(sample_name, "_", CHROM, ".CpG_report.txt.gz"))
}

# Pool a list of replicate sample names using DMRcaller's readBismarkPool,
# which sums readsM and readsN across replicates per CpG site (identical to
# what the original build_chr1_cache.R did).
pool_reps <- function(sample_names) {
  paths <- vapply(sample_names, rep_path, character(1))
  stopifnot(all(file.exists(paths)))
  readBismarkPool(paths)
}

# Enumerate the 10 unique 3-vs-3 partitions of 6 replicates.
# Returns a data.frame with one row per partition:
#   perm_id  : 0..9 (0 = real / identity)
#   groupA   : comma-separated sample names assigned to group A
#   groupB   : comma-separated sample names assigned to group B
#   is_real  : TRUE only for perm_id == 0
enumerate_partitions <- function() {
  idx_A <- combn(6, 3)               # 6 x 20 matrix of group-A indices
  partitions <- list()
  seen_keys  <- character(0)
  for (k in seq_len(ncol(idx_A))) {
    a <- sort(idx_A[, k]);  b <- sort(setdiff(1:6, a))
    # canonical key: smaller side first to dedupe mirror partitions
    if (a[1] < b[1]) { side1 <- a; side2 <- b } else { side1 <- b; side2 <- a }
    key <- paste(c(side1, side2), collapse = "-")
    if (key %in% seen_keys) next
    seen_keys <- c(seen_keys, key)
    partitions[[length(partitions) + 1]] <- list(A = a, B = b)
  }
  # 10 unique partitions expected
  stopifnot(length(partitions) == 10)
  # Identify the identity (real) partition: group A == ASO_VPA reps (1:3)
  is_identity <- vapply(partitions, function(p) all(sort(p$A) == 1:3), logical(1))
  stopifnot(sum(is_identity) == 1)
  # Reorder so identity is perm_id = 0 (first row)
  ord <- c(which(is_identity), which(!is_identity))
  partitions <- partitions[ord]
  out <- data.frame(
    perm_id = 0:9,
    groupA  = vapply(partitions, function(p) paste(ALL_6[p$A], collapse = ","), character(1)),
    groupB  = vapply(partitions, function(p) paste(ALL_6[p$B], collapse = ","), character(1)),
    is_real = c(TRUE, rep(FALSE, 9)),
    stringsAsFactors = FALSE
  )
  out
}

# Pre-build the pooled GRanges for one partition. Returns list(treat=, ctrl=).
build_pair_for_perm <- function(perm_row) {
  a_samples <- strsplit(perm_row$groupA, ",", fixed = TRUE)[[1]]
  b_samples <- strsplit(perm_row$groupB, ",", fixed = TRUE)[[1]]
  list(
    treat = pool_reps(a_samples),
    ctrl  = pool_reps(b_samples)
  )
}

# threshold lookup (matches parameter_benchmark_helpers.R::get_thresholds, strict=TRUE)
get_thresholds_strict <- function(method, ws = 300) {
  list(
    pval    = 0.01,
    minCpG  = 4,
    minDiff = if (method == "noise_filter") 0.4 else 0.2,
    minSize = if (method == "neighbourhood") 1 else 50,
    minGap  = if (method == "noise_filter") 0 else 200,
    test    = if (method == "noise_filter") "score" else "fisher"
  )
}

# single DMR call (matches run_dmrs_one logic)
run_dmr_call <- function(treat, ctrl, method, ws, kernel = "triangular", region) {
  th <- get_thresholds_strict(method, ws)
  setTimeLimit(cpu = 7200, elapsed = 7200, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE))
  tryCatch({
    if (method == "noise_filter") {
      computeDMRs(treat, ctrl, regions = region, context = "CG",
        method = "noise_filter", windowSize = ws, kernelFunction = kernel,
        test = th$test, pValueThreshold = th$pval, minCytosinesCount = th$minCpG,
        minProportionDifference = th$minDiff, minGap = th$minGap, minSize = th$minSize,
        minReadsPerCytosine = 4, parallel = FALSE)
    } else if (method == "bins") {
      computeDMRs(treat, ctrl, regions = region, context = "CG",
        method = "bins", binSize = ws, test = th$test, pValueThreshold = th$pval,
        minCytosinesCount = th$minCpG, minProportionDifference = th$minDiff,
        minGap = th$minGap, minSize = th$minSize, minReadsPerCytosine = 4,
        parallel = FALSE)
    } else {
      computeDMRs(treat, ctrl, regions = region, context = "CG",
        method = "neighbourhood", test = th$test, pValueThreshold = th$pval,
        minCytosinesCount = th$minCpG, minProportionDifference = th$minDiff,
        minGap = th$minGap, minSize = th$minSize, minReadsPerCytosine = 4,
        parallel = FALSE)
    }
  }, error = function(e) {
    message("  Error: ", e$message)
    structure(list(), class = "dmr_error", message = e$message)
  })
}

is_error <- function(x) inherits(x, "dmr_error")

dmr_size_summary <- function(dmrs) {
  if (is_error(dmrs) || length(dmrs) == 0)
    return(list(median = NA_real_, mean = NA_real_, q1 = NA_real_,
                q3 = NA_real_, min = NA_real_, max = NA_real_))
  w <- as.numeric(width(dmrs))
  list(median = median(w), mean = round(mean(w), 1),
       q1 = as.numeric(quantile(w, 0.25)),
       q3 = as.numeric(quantile(w, 0.75)),
       min = min(w), max = max(w))
}
