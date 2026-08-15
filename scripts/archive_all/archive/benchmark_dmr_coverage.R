.libPaths("~/R/library")
# benchmark_dmr_coverage.R
# Runs bins and noise_filter on real chr1 data only
# Extracts exact DMR widths (sum of bp covered) per parameter combination
# Used to generate accurate genome coverage plots (Radu-style Panel A/B)
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
REGION_END <- 248956422

message("Loading chr1 data...")
aso_vpa <- readBismarkPool(c(
  file.path(COV_DIR, "ASO_VPA_1_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_2_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_3_chr1.CpG_report.txt.gz")
))
aso_ctrl <- readBismarkPool(c(
  file.path(COV_DIR, "ASO_CTRL_1_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_CTRL_2_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_CTRL_3_chr1.CpG_report.txt.gz")
))
message("Loaded. ASO_VPA CpGs: ", length(aso_vpa))

window_sizes <- c(100, 200, 300, 500, 1000, 2000)
results <- list()

get_coverage <- function(dmrs) {
  if (is.null(dmrs) || length(dmrs) == 0) return(list(n=0, bp=0))
  list(n=length(dmrs), bp=sum(width(dmrs)))
}

run_method <- function(treat, ctrl, method, ws, kernel="triangular", strict=FALSE) {
  th <- if (strict) {
    list(pval=0.01, minCpG=4, minDiff=0.2, minSize=50)
  } else {
    list(pval=0.05, minCpG=1, minDiff=0.1, minSize=1)
  }
  regions <- GRanges("chr1", IRanges(seq(1, REGION_END, by=ws), width=ws))
  tryCatch({
    if (method == "bins") {
      computeDMRs(treat, ctrl, regions=regions,
        context="CG", method="bins", binSize=ws,
        test="fisher", pValueThreshold=th$pval,
        minCytosinesCount=th$minCpG,
        minProportionDifference=th$minDiff,
        minGap=0, minSize=th$minSize,
        minReadsPerCytosine=4, cores=36, parallel=TRUE)
    } else {
      computeDMRs(treat, ctrl, regions=regions,
        context="CG", method="noise_filter",
        windowSize=ws, kernelFunction=kernel,
        test="score", pValueThreshold=th$pval,
        minCytosinesCount=th$minCpG,
        minProportionDifference=th$minDiff,
        minGap=0, minSize=th$minSize,
        minReadsPerCytosine=4, cores=36, parallel=TRUE)
    }
  }, error=function(e) { message("Error: ", e$message); NULL })
}

for (ws in window_sizes) {
  for (strict in c(FALSE, TRUE)) {
    mode <- if (strict) "strict" else "loose"

    message("\n--- bins | ws=", ws, " | ", mode, " ---")
    dmrs <- run_method(aso_vpa, aso_ctrl, "bins", ws, strict=strict)
    cov  <- get_coverage(dmrs)
    message("  n=", cov$n, " bp=", cov$bp)
    results[[paste("bins", ws, mode)]] <- data.frame(
      method="bins", window_size=ws, mode=mode,
      kernel="NA", n_dmrs=cov$n, total_bp=cov$bp)

    message("--- noise_filter | ws=", ws, " | ", mode, " ---")
    dmrs <- run_method(aso_vpa, aso_ctrl, "noise_filter", ws,
                       kernel="triangular", strict=strict)
    cov  <- get_coverage(dmrs)
    message("  n=", cov$n, " bp=", cov$bp)
    results[[paste("nf", ws, mode)]] <- data.frame(
      method="noise_filter", window_size=ws, mode=mode,
      kernel="triangular", n_dmrs=cov$n, total_bp=cov$bp)
  }
}

df <- do.call(rbind, results)
write.csv(df, file.path(OUT_DIR, "benchmark_dmr_coverage_real.csv"), row.names=FALSE)
message("\n=== COVERAGE RESULTS ===")
print(df, row.names=FALSE)
message("\nSaved to: ", OUT_DIR, "/benchmark_dmr_coverage_real.csv")
