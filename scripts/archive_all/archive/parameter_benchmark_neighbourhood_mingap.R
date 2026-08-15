.libPaths("~/R/library")
# parameter_benchmark_neighbourhood_mingap.R
# Neighbourhood minGap benchmark — label-swap null model
# minGap varied at 100, 200, 300, 500, 1000, 2000bp on whole chr1
# Complements main benchmark — minGap is the relevant parameter for neighbourhood
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(ggplot2)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 248956422
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chr1_region <- GRanges("chr1", IRanges(1, REGION_END))

load_group_chr1 <- function(samples) {
  paths <- file.path(COV_DIR, paste0(samples, "_", CHROM, ".CpG_report.txt.gz"))
  paths <- paths[file.exists(paths)]
  message("  Loading ", length(paths), " files...")
  readBismarkPool(paths)
}

message("Loading chr1 data...")
aso_vpa  <- load_group_chr1(c("ASO_VPA_1", "ASO_VPA_2", "ASO_VPA_3"))
aso_ctrl <- load_group_chr1(c("ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3"))
message("  ASO_VPA CpGs:  ", length(aso_vpa))
message("  ASO_CTRL CpGs: ", length(aso_ctrl))

# Label-swap null model
message("Creating label-swap scrambled datasets...")
aso_vpa_scr  <- aso_ctrl
aso_ctrl_scr <- aso_vpa
message("  Done")

get_thresholds <- function(strict) {
  if (!strict) {
    list(pval=0.05, minCpG=1, minDiff=0.1, minSize=1, test="fisher")
  } else {
    list(pval=0.01, minCpG=4, minDiff=0.2, minSize=50, test="fisher")
  }
}

run_neighbourhood <- function(treat, ctrl, mingap, strict=FALSE) {
  th <- get_thresholds(strict)
  tryCatch({
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
      cores                   = 32,
      parallel                = TRUE)
  }, error = function(e) { message("  Error: ", e$message); NULL })
}

mingap_sizes <- c(100, 200, 300, 500, 1000, 2000)
modes        <- c(FALSE, TRUE)
results      <- list()

for (mg in mingap_sizes) {
  for (strict in modes) {
    mode_lab <- if (strict) "strict" else "loose"
    message("\n--- neighbourhood_mingap | minGap=", mg, " | ", mode_lab, " ---")

    message("  Real data...")
    dmrs_real <- run_neighbourhood(aso_vpa, aso_ctrl, mg, strict)
    n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
    message("  Real DMRs: ", n_real)

    message("  Scrambled...")
    dmrs_scr <- run_neighbourhood(aso_vpa_scr, aso_ctrl_scr, mg, strict)
    n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
    message("  Scrambled DMRs: ", n_scr)

    ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) {
      round(n_real / n_scr, 2)
    } else if (!is.na(n_real) && !is.na(n_scr) && n_scr == 0) {
      Inf
    } else NA
    message("  Signal/Noise: ", ratio)

    results[[paste("neighbourhood_mingap", mg, mode_lab, sep="_")]] <- data.frame(
      method          = "neighbourhood_mingap",
      window_size     = mg,
      mode            = mode_lab,
      kernel          = "NA",
      scramble_method = "label_swap",
      n_real          = n_real,
      n_scrambled     = n_scr,
      ratio           = ratio,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, results)
write.csv(summary_df,
          file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap_label_swap.csv"),
          row.names = FALSE)

message("\n=== NEIGHBOURHOOD minGap RESULTS (Label-swap) ===")
message(sprintf("%-8s %-6s %-12s %-12s %-10s",
                "Mode", "minGap", "Real DMRs", "Scr DMRs", "S/N"))
message(paste(rep("-", 50), collapse=""))
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i,]
  message(sprintf("%-8s %-6s %-12s %-12s %-10s",
                  r$mode, r$window_size, r$n_real, r$n_scrambled, r$ratio))
}
message("\nDone. Saved to: ", OUT_DIR,
        "/parameter_benchmark_neighbourhood_mingap_label_swap.csv")
