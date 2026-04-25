.libPaths("~/R/library")
# parameter_benchmark_noisefilter.R
# noise_filter benchmark — 100, 200, 300, 500bp on chr1:1-10Mb
# Run separately due to compute time — submit with 72h, 128G, 20 cores
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(ggplot2)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 10000000

load_group_chr1 <- function(samples) {
  paths <- file.path(COV_DIR, paste0(samples, "_", CHROM, ".CpG_report.txt.gz"))
  paths <- paths[file.exists(paths)]
  message("  Loading ", length(paths), " files...")
  dat <- readBismarkPool(paths)
  dat[start(dat) <= REGION_END]
}

message("Loading chr1 data (first 10Mb)...")
aso_vpa  <- load_group_chr1(c("ASO_VPA_1",  "ASO_VPA_2",  "ASO_VPA_3"))
aso_ctrl <- load_group_chr1(c("ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3"))

scramble_counts <- function(dat) {
  set.seed(42)
  dat_scr <- dat
  idx <- sample(length(dat))
  mcols(dat_scr)$readsM <- mcols(dat)$readsM[idx]
  mcols(dat_scr)$readsN <- mcols(dat)$readsN[idx]
  dat_scr
}

message("Generating scrambled datasets...")
aso_vpa_scr  <- scramble_counts(aso_vpa)
aso_ctrl_scr <- scramble_counts(aso_ctrl)

make_tiles <- function(binsize) {
  GRanges("chr1", IRanges(
    start = seq(1, REGION_END, by = binsize),
    width = binsize))
}

run_noisefilter <- function(treat, ctrl, ws, kernel, strict = FALSE) {
  regions <- make_tiles(ws)

  if (!strict) {
    pval <- 0.05; minCpG <- 1; minDiff <- 0.1; minSize <- 1
    test <- "fisher"
  } else {
    pval <- 0.01; minCpG <- 4; minDiff <- 0.4; minSize <- 50
    test <- "score"
  }

  tryCatch({
    computeDMRs(treat, ctrl,
      regions                 = regions,
      context                 = "CG",
      method                  = "noise_filter",
      windowSize              = ws,
      kernelFunction          = kernel,
      test                    = test,
      pValueThreshold         = pval,
      minCytosinesCount       = minCpG,
      minProportionDifference = minDiff,
      minGap                  = 0,
      minSize                 = minSize,
      minReadsPerCytosine     = 4,
      cores                   = 20)
  }, error = function(e) { message("  Error: ", e$message); NULL })
}

window_sizes <- c(100, 200, 300, 500, 1000, 2000)
kernels      <- c("triangular", "uniform", "epanechnicov")
modes        <- c(FALSE, TRUE)
results      <- list()

for (ws in window_sizes) {
  for (ker in kernels) {
    for (strict in modes) {
      mode_lab <- if (strict) "strict" else "loose"
      message("\n--- noise_filter | ws=", ws, " | kernel=", ker, " | ", mode_lab, " ---")

      message("  Real data...")
      dmrs_real <- run_noisefilter(aso_vpa, aso_ctrl, ws, ker, strict)
      n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
      message("  Real DMRs: ", n_real)

      message("  Scrambled...")
      dmrs_scr <- run_noisefilter(aso_vpa_scr, aso_ctrl_scr, ws, ker, strict)
      n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
      message("  Scrambled DMRs: ", n_scr)

      ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) {
        round(n_real / n_scr, 2)
      } else if (!is.na(n_real) && !is.na(n_scr) && n_scr == 0) {
        Inf
      } else NA
      message("  Signal/Noise: ", ratio)

      results[[paste("noise_filter", ws, ker, mode_lab, sep="_")]] <- data.frame(
        method      = "noise_filter",
        window_size = ws,
        mode        = mode_lab,
        kernel      = ker,
        n_real      = n_real,
        n_scrambled = n_scr,
        ratio       = ratio,
        stringsAsFactors = FALSE
      )
    }
  }
}

noisefilter_df <- do.call(rbind, results)
write.csv(noisefilter_df,
          file.path(OUT_DIR, "parameter_benchmark_noisefilter.csv"),
          row.names = FALSE)
message("\nDone. Results saved to: ", OUT_DIR, "/parameter_benchmark_noisefilter.csv")
