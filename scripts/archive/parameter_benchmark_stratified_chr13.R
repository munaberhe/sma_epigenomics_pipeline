.libPaths("~/R/library")
# parameter_benchmark_stratified.R
# DMR parameter benchmarking — stratified scramble null model
# Scramble preserves coverage distribution by shuffling within coverage strata
# and uses same permutation index for both datasets to avoid systematic differences
# Identical structure to label-swap and Archie scripts
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(ggplot2)

COV_DIR  <- "results/alignments/bs/by_chr"
OUT_DIR  <- "results/dmr_benchmark"
CHROM    <- "chr13"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REGION_END <- 114364328

load_group_chr1 <- function(samples) {
  paths <- file.path(COV_DIR, paste0(samples, "_", CHROM, ".CpG_report.txt.gz"))
  paths <- paths[file.exists(paths)]
  message("  Loading ", length(paths), " files...")
  readBismarkPool(paths)
}

message("Loading chr1 data...")
aso_vpa  <- load_group_chr1(c("ASO_VPA_1",  "ASO_VPA_2",  "ASO_VPA_3"))
aso_ctrl <- load_group_chr1(c("ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3"))
message("  ASO_VPA CpGs:  ", length(aso_vpa))
message("  ASO_CTRL CpGs: ", length(aso_ctrl))

# Stratified scramble — shuffle within coverage strata using same index for both
# Fixes: (1) coverage distribution mismatch, (2) independent seed artefact
scramble_stratified <- function(dat1, dat2, seed=42) {
  set.seed(seed)
  coverage <- mcols(dat1)$readsN
  strata <- cut(coverage,
                breaks = c(0, 5, 10, 20, 50, Inf),
                labels = c("1-5", "6-10", "11-20", "21-50", "50+"),
                include.lowest = TRUE)

  idx <- seq_along(coverage)
  for (s in levels(strata)) {
    stratum_idx <- which(strata == s)
    if (length(stratum_idx) > 1) {
      idx[stratum_idx] <- sample(stratum_idx)
    }
  }

  dat1_scr <- dat1
  dat2_scr <- dat2
  mcols(dat1_scr)$readsM <- mcols(dat1)$readsM[idx]
  mcols(dat1_scr)$readsN <- mcols(dat1)$readsN[idx]
  mcols(dat2_scr)$readsM <- mcols(dat2)$readsM[idx]
  mcols(dat2_scr)$readsN <- mcols(dat2)$readsN[idx]

  message("  Strata sizes: ", paste(table(strata), collapse=", "))
  message("  Original ASO_VPA mean methylation:  ",
          round(sum(mcols(dat1)$readsM)/sum(mcols(dat1)$readsN)*100, 1), "%")
  message("  Scrambled ASO_VPA mean methylation: ",
          round(sum(mcols(dat1_scr)$readsM)/sum(mcols(dat1_scr)$readsN)*100, 1), "%")

  list(dat1=dat1_scr, dat2=dat2_scr)
}

message("Creating stratified scrambled datasets...")
scrambled <- scramble_stratified(aso_vpa, aso_ctrl, seed=42)
aso_vpa_scr  <- scrambled$dat1
aso_ctrl_scr <- scrambled$dat2
message("  Done")

chr_region <- GRanges("chr13", IRanges(1, REGION_END))

make_tiles <- function(binsize) {
  GRanges("chr13", IRanges(
    start = seq(1, REGION_END, by = binsize),
    width = binsize))
}

get_thresholds <- function(method, strict) {
  if (!strict) {
    list(pval=0.05, minCpG=1, minDiff=0.1, minSize=1, test="fisher")
  } else {
    list(pval=0.01, minCpG=4,
         minDiff=if (method == "noise_filter") 0.4 else 0.2,
         minSize=50,
         test=if (method == "noise_filter") "score" else "fisher")
  }
}

run_dmrs <- function(treat, ctrl, method, ws, kernel="triangular", strict=FALSE) {
  th      <- get_thresholds(method, strict)
  regions <- make_tiles(ws)

  tryCatch({
    if (method == "noise_filter") {
      computeDMRs(treat, ctrl,
        regions                 = regions,
        context                 = "CG",
        method                  = "noise_filter",
        windowSize              = ws,
        kernelFunction          = kernel,
        test                    = th$test,
        pValueThreshold         = th$pval,
        minCytosinesCount       = th$minCpG,
        minProportionDifference = th$minDiff,
        minGap                  = 0,
        minSize                 = th$minSize,
        minReadsPerCytosine     = 4,
        cores                   = 32,
        parallel                = TRUE)
    } else if (method == "bins") {
      computeDMRs(treat, ctrl,
        regions                 = regions,
        context                 = "CG",
        method                  = "bins",
        binSize                 = ws,
        test                    = th$test,
        pValueThreshold         = th$pval,
        minCytosinesCount       = th$minCpG,
        minProportionDifference = th$minDiff,
        minGap                  = 0,
        minSize                 = th$minSize,
        minReadsPerCytosine     = 4,
        cores                   = 32,
        parallel                = TRUE)
    } else {
      computeDMRs(treat, ctrl,
        regions                 = regions,
        context                 = "CG",
        method                  = "neighbourhood",
        test                    = th$test,
        pValueThreshold         = th$pval,
        minCytosinesCount       = th$minCpG,
        minProportionDifference = th$minDiff,
        minGap                  = 0,
        minSize                 = th$minSize,
        minReadsPerCytosine     = 4,
        cores                   = 32,
        parallel                = TRUE)
    }
  }, error = function(e) { message("  Error: ", e$message); NULL })
}

window_sizes <- c(100, 200, 300, 500, 1000, 2000)
methods      <- c("bins", "neighbourhood", "noise_filter")
kernels_nf   <- c("triangular", "uniform", "epanechnicov")
modes        <- c(FALSE, TRUE)
results      <- list()

for (method in methods) {
  for (ws in window_sizes) {
    ker_vec <- if (method == "noise_filter") kernels_nf else "NA"
    for (strict in modes) {
      for (ker in ker_vec) {
        mode_lab <- if (strict) "strict" else "loose"
        message("\n--- ", mode_lab, " | ", method, " | ws=", ws, " | kernel=", ker, " ---")

        message("  Real data...")
        dmrs_real <- run_dmrs(aso_vpa, aso_ctrl, method, ws, kernel=ker, strict=strict)
        n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
        message("  Real DMRs: ", n_real)

        message("  Scrambled...")
        dmrs_scr <- run_dmrs(aso_vpa_scr, aso_ctrl_scr, method, ws, kernel=ker, strict=strict)
        n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
        message("  Scrambled DMRs: ", n_scr)

        ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) {
          round(n_real / n_scr, 2)
        } else if (!is.na(n_real) && !is.na(n_scr) && n_scr == 0) {
          Inf
        } else NA
        message("  Signal/Noise: ", ratio)

        results[[paste(method, ws, mode_lab, ker, sep="_")]] <- data.frame(
          method          = method,
          window_size     = ws,
          mode            = mode_lab,
          kernel          = ker,
          scramble_method = "stratified_scramble",
          n_real          = n_real,
          n_scrambled     = n_scr,
          ratio           = ratio,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

summary_df <- do.call(rbind, results)
write.csv(summary_df,
          file.path(OUT_DIR, "parameter_benchmark_stratified_scramble_chr13.csv"),
          row.names = FALSE)

message("\n=== BENCHMARK RESULTS (Stratified scramble) ===")
message(sprintf("%-8s %-15s %-6s %-14s %-12s %-12s %-10s",
                "Mode", "Method", "ws", "Kernel", "Real DMRs", "Scr DMRs", "S/N"))
message(paste(rep("-", 80), collapse=""))
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i,]
  message(sprintf("%-8s %-15s %-6s %-14s %-12s %-12s %-10s",
                  r$mode, r$method, r$window_size, r$kernel,
                  r$n_real, r$n_scrambled, r$ratio))
}

message("\nDone. Saved to: ", OUT_DIR,
        "/parameter_benchmark_stratified_scramble_chr13.csv")
