.libPaths("~/R/library")
# parameter_benchmark_archie.R
# DMR parameter benchmarking — Archie scramble null model
# All methods (bins, neighbourhood, noise_filter) tested at all window sizes
# Two threshold regimes (loose/strict) and three kernels for noise_filter
# Identical structure to parameter_benchmark.R — only scrambling method differs
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(ggplot2)

COV_DIR  <- "results/alignments/bs/by_chr"
OUT_DIR  <- "results/dmr_benchmark"
CHROM    <- "chr1"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REGION_END <- 248956422  # full chr1

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

# Archie scramble null model — shuffle readsM/readsN across positions
# breaks spatial signal while keeping coverage distribution intact
scramble_counts <- function(dat, seed=42) {
  set.seed(seed)
  dat_scr <- dat
  idx <- sample(length(dat))
  mcols(dat_scr)$readsM <- mcols(dat)$readsM[idx]
  mcols(dat_scr)$readsN <- mcols(dat)$readsN[idx]
  dat_scr
}

message("Creating Archie scrambled datasets...")
aso_vpa_scr  <- scramble_counts(aso_vpa,  seed=42)
aso_ctrl_scr <- scramble_counts(aso_ctrl, seed=123)
message("  Done")

chr1_region <- GRanges("chr1", IRanges(1, REGION_END))

make_tiles <- function(binsize) {
  GRanges("chr1", IRanges(
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
      # neighbourhood — tiled so window size actually varies results
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
          scramble_method = "archie_scramble",
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
          file.path(OUT_DIR, "parameter_benchmark_archie_scramble.csv"),
          row.names = FALSE)

message("\n=== BENCHMARK RESULTS (Archie scramble) ===")
message(sprintf("%-8s %-15s %-6s %-14s %-12s %-12s %-10s",
                "Mode", "Method", "ws", "Kernel", "Real DMRs", "Scr DMRs", "S/N"))
message(paste(rep("-", 80), collapse=""))
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i,]
  message(sprintf("%-8s %-15s %-6s %-14s %-12s %-12s %-10s",
                  r$mode, r$method, r$window_size, r$kernel,
                  r$n_real, r$n_scrambled, r$ratio))
}

message("\nDone. Results saved to: ", OUT_DIR,
        "/parameter_benchmark_archie_scramble.csv")
