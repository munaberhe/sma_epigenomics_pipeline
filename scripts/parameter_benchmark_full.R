.libPaths("~/R/library")
# parameter_benchmark_full.R
# Full DMR parameter benchmark — bins, neighbourhood (windowSize + minGap), noise_filter
# All methods tested at 100, 200, 300, 500, 1000, 2000bp on chr1:1-10Mb
# Scrambling based on Archie's approach — shuffles readsM/readsN to break spatial signal
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(ggplot2)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 10000000
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chr1_region <- GRanges("chr1", IRanges(1, REGION_END))

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
message("  ASO_VPA CpGs:  ", length(aso_vpa))
message("  ASO_CTRL CpGs: ", length(aso_ctrl))

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
message("  Done")

make_tiles <- function(binsize) {
  GRanges("chr1", IRanges(
    start = seq(1, REGION_END, by = binsize),
    width = binsize))
}

# shared threshold logic
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

run_dmrs <- function(treat, ctrl, method, ws, kernel="triangular",
                     strict=FALSE, mingap=0) {
  th <- get_thresholds(method, strict)

  tryCatch({
    if (method == "bins") {
      computeDMRs(treat, ctrl,
        regions                 = make_tiles(ws),
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
        cores                   = 20)

    } else if (method == "neighbourhood") {
      computeDMRs(treat, ctrl,
        regions                 = make_tiles(ws),
        context                 = "CG",
        method                  = "neighbourhood",
        test                    = th$test,
        pValueThreshold         = th$pval,
        minCytosinesCount       = th$minCpG,
        minProportionDifference = th$minDiff,
        minGap                  = 0,
        minSize                 = th$minSize,
        minReadsPerCytosine     = 4,
        cores                   = 20)

    } else if (method == "neighbourhood_mingap") {
      computeDMRs(treat, ctrl,
        regions                 = chr1_region,
        context                 = "CG",
        method                  = "neighbourhood",
        test                    = th$test,
        pValueThreshold         = th$pval,
        minCytosinesCount       = th$minCpG,
        minProportionDifference = th$minDiff,
        minGap                  = mingap,
        minSize                 = th$minSize,
        minReadsPerCytosine     = 4,
        cores                   = 20)

    } else if (method == "noise_filter") {
      computeDMRs(treat, ctrl,
        regions                 = make_tiles(ws),
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
        cores                   = 20)
    }
  }, error = function(e) { message("  Error: ", e$message); NULL })
}

record_result <- function(results, key, method, ws, mode, kernel, n_real, n_scr) {
  ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) {
    round(n_real / n_scr, 2)
  } else if (!is.na(n_real) && !is.na(n_scr) && n_scr == 0) {
    Inf
  } else NA
  message("  Signal/Noise: ", ratio)
  results[[key]] <- data.frame(
    method=method, window_size=ws, mode=mode, kernel=kernel,
    n_real=n_real, n_scrambled=n_scr, ratio=ratio,
    stringsAsFactors=FALSE)
  results
}

sizes   <- c(100, 200, 300, 500, 1000, 2000)
modes   <- c(FALSE, TRUE)
kernels <- c("triangular", "uniform", "epanechnicov")
results <- list()

# --- BINS ---
message("\n====== BINS ======")
for (ws in sizes) {
  for (strict in modes) {
    mode_lab <- if (strict) "strict" else "loose"
    message("\n--- bins | ws=", ws, " | ", mode_lab, " ---")
    dmrs_real <- run_dmrs(aso_vpa, aso_ctrl, "bins", ws, strict=strict)
    n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
    message("  Real DMRs: ", n_real)
    dmrs_scr <- run_dmrs(aso_vpa_scr, aso_ctrl_scr, "bins", ws, strict=strict)
    n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
    message("  Scrambled DMRs: ", n_scr)
    results <- record_result(results, paste("bins", ws, mode_lab, sep="_"),
                             "bins", ws, mode_lab, "NA", n_real, n_scr)
  }
}

# --- NEIGHBOURHOOD (windowSize) ---
message("\n====== NEIGHBOURHOOD (windowSize) ======")
for (ws in sizes) {
  for (strict in modes) {
    mode_lab <- if (strict) "strict" else "loose"
    message("\n--- neighbourhood | ws=", ws, " | ", mode_lab, " ---")
    dmrs_real <- run_dmrs(aso_vpa, aso_ctrl, "neighbourhood", ws, strict=strict)
    n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
    message("  Real DMRs: ", n_real)
    dmrs_scr <- run_dmrs(aso_vpa_scr, aso_ctrl_scr, "neighbourhood", ws, strict=strict)
    n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
    message("  Scrambled DMRs: ", n_scr)
    results <- record_result(results, paste("neighbourhood", ws, mode_lab, sep="_"),
                             "neighbourhood", ws, mode_lab, "NA", n_real, n_scr)
  }
}

# --- NEIGHBOURHOOD (minGap) ---
message("\n====== NEIGHBOURHOOD (minGap) ======")
for (mg in sizes) {
  for (strict in modes) {
    mode_lab <- if (strict) "strict" else "loose"
    message("\n--- neighbourhood_mingap | minGap=", mg, " | ", mode_lab, " ---")
    dmrs_real <- run_dmrs(aso_vpa, aso_ctrl, "neighbourhood_mingap",
                          ws=100, strict=strict, mingap=mg)
    n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
    message("  Real DMRs: ", n_real)
    dmrs_scr <- run_dmrs(aso_vpa_scr, aso_ctrl_scr, "neighbourhood_mingap",
                         ws=100, strict=strict, mingap=mg)
    n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
    message("  Scrambled DMRs: ", n_scr)
    results <- record_result(results, paste("neighbourhood_mingap", mg, mode_lab, sep="_"),
                             "neighbourhood_mingap", mg, mode_lab, "NA", n_real, n_scr)
  }
}

# --- NOISE_FILTER ---
message("\n====== NOISE_FILTER ======")
for (ws in sizes) {
  for (ker in kernels) {
    for (strict in modes) {
      mode_lab <- if (strict) "strict" else "loose"
      message("\n--- noise_filter | ws=", ws, " | kernel=", ker, " | ", mode_lab, " ---")
      dmrs_real <- run_dmrs(aso_vpa, aso_ctrl, "noise_filter",
                            ws, kernel=ker, strict=strict)
      n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
      message("  Real DMRs: ", n_real)
      dmrs_scr <- run_dmrs(aso_vpa_scr, aso_ctrl_scr, "noise_filter",
                           ws, kernel=ker, strict=strict)
      n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
      message("  Scrambled DMRs: ", n_scr)
      results <- record_result(results,
                               paste("noise_filter", ws, ker, mode_lab, sep="_"),
                               "noise_filter", ws, mode_lab, ker, n_real, n_scr)
    }
  }
}

summary_df <- do.call(rbind, results)
write.csv(summary_df,
          file.path(OUT_DIR, "parameter_benchmark_full.csv"),
          row.names=FALSE)

message("\n=== BENCHMARK COMPLETE ===")
message("Saved to: ", OUT_DIR, "/parameter_benchmark_full.csv")
message("Total rows: ", nrow(summary_df))
