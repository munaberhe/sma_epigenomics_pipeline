.libPaths("~/R/library")
# parameter_benchmark_neighbourhood_mingap.R
# Neighbourhood method benchmark — varying minGap (0, 100, 200, 500bp) on chr1:1-10Mb
# minGap is the relevant parameter for neighbourhood (not windowSize)
# Saves to separate file — does not overwrite parameter_benchmark_neighbourhood.csv
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(ggplot2)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 10000000
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

run_neighbourhood_mingap <- function(treat, ctrl, mingap, strict = FALSE) {
  if (!strict) {
    pval <- 0.05; minCpG <- 1; minDiff <- 0.1; minSize <- 1
  } else {
    pval <- 0.01; minCpG <- 4; minDiff <- 0.2; minSize <- 50
  }

  tryCatch({
    computeDMRs(treat, ctrl,
      regions                 = chr1_region,
      context                 = "CG",
      method                  = "neighbourhood",
      test                    = "fisher",
      pValueThreshold         = pval,
      minCytosinesCount       = minCpG,
      minProportionDifference = minDiff,
      minGap                  = mingap,
      minSize                 = minSize,
      minReadsPerCytosine     = 4,
      cores                   = 20)
  }, error = function(e) { message("  Error: ", e$message); NULL })
}

# varying minGap instead of windowSize — this is the correct parameter for neighbourhood
mingap_sizes <- c(0, 100, 200, 500, 1000, 2000)
modes        <- c(FALSE, TRUE)
results      <- list()

for (mg in mingap_sizes) {
  for (strict in modes) {
    mode_lab <- if (strict) "strict" else "loose"
    message("\n--- neighbourhood_mingap | minGap=", mg, " | ", mode_lab, " ---")

    message("  Real data...")
    dmrs_real <- run_neighbourhood_mingap(aso_vpa, aso_ctrl, mg, strict)
    n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
    message("  Real DMRs: ", n_real)

    message("  Scrambled...")
    dmrs_scr <- run_neighbourhood_mingap(aso_vpa_scr, aso_ctrl_scr, mg, strict)
    n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
    message("  Scrambled DMRs: ", n_scr)

    ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) {
      round(n_real / n_scr, 2)
    } else if (!is.na(n_real) && !is.na(n_scr) && n_scr == 0) {
      Inf
    } else NA
    message("  Signal/Noise: ", ratio)

    results[[paste("neighbourhood_mingap", mg, mode_lab, sep="_")]] <- data.frame(
      method      = "neighbourhood_mingap",
      window_size = mg,
      mode        = mode_lab,
      kernel      = "NA",
      n_real      = n_real,
      n_scrambled = n_scr,
      ratio       = ratio,
      stringsAsFactors = FALSE
    )
  }
}

neighbourhood_mingap_df <- do.call(rbind, results)
write.csv(neighbourhood_mingap_df,
          file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap.csv"),
          row.names = FALSE)
message("\nDone. Results saved to: ", OUT_DIR, "/parameter_benchmark_neighbourhood_mingap.csv")
