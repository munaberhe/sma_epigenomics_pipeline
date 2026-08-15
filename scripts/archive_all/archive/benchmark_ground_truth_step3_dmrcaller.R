.libPaths("~/R/library")
# benchmark_ground_truth_step3_dmrcaller.R
# Step 3: Run all DMR calling methods on spiked datasets
# Tests bins, neighbourhood, noise_filter at multiple window sizes
# against ground truth spike-in regions
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

OUT_DIR   <- "results/benchmark_ground_truth"
SPIKE_CSV <- file.path(OUT_DIR, "spike_in_regions.csv")

# Load spike-in regions — these are the ground truth
spike_df <- read.csv(SPIKE_CSV)
spike_regions <- GRanges(
  seqnames = spike_df$chr,
  ranges   = IRanges(spike_df$start, spike_df$end)
)
message("Ground truth: ", length(spike_regions), " spike-in regions")

# Load original control — never modified
message("Loading control data...")
aso_ctrl <- readRDS(file.path(OUT_DIR, "aso_ctrl_original.rds"))

# Parameters to test
window_sizes <- c(100, 200, 300, 500, 1000, 2000)
delta_meths  <- c(0.1, 0.2, 0.3, 0.4)
methods      <- c("bins", "neighbourhood", "noise_filter")
kernels_nf   <- c("triangular")
chr1_region  <- GRanges("chr1", IRanges(1, 248956422))

# Function to run DMRcaller and evaluate against ground truth
run_and_evaluate <- function(treatment, control, method, ws,
                             kernel="triangular", strict=FALSE) {
  th <- if (!strict) {
    list(pval=0.05, minCpG=1, minDiff=0.1, minSize=1, test="fisher")
  } else {
    list(pval=0.01, minCpG=4, minDiff=0.2, minSize=50, test="fisher")
  }

  dmrs <- tryCatch({
    if (method == "bins") {
      computeDMRs(treatment, control,
        regions=GRanges("chr1", IRanges(seq(1,248956422,by=ws), width=ws)),
        context="CG", method="bins", binSize=ws,
        test=th$test, pValueThreshold=th$pval,
        minCytosinesCount=th$minCpG,
        minProportionDifference=th$minDiff,
        minGap=0, minSize=th$minSize,
        minReadsPerCytosine=4, cores=64, parallel=TRUE)
    } else if (method == "neighbourhood") {
      computeDMRs(treatment, control,
        regions=chr1_region,
        context="CG", method="neighbourhood",
        test=th$test, pValueThreshold=th$pval,
        minCytosinesCount=th$minCpG,
        minProportionDifference=th$minDiff,
        minGap=0, minSize=th$minSize,
        minReadsPerCytosine=4, cores=64, parallel=TRUE)
    } else {
      computeDMRs(treatment, control,
        regions=GRanges("chr1", IRanges(seq(1,248956422,by=ws), width=ws)),
        context="CG", method="noise_filter",
        windowSize=ws, kernelFunction=kernel,
        test=th$test, pValueThreshold=th$pval,
        minCytosinesCount=th$minCpG,
        minProportionDifference=th$minDiff,
        minGap=0, minSize=th$minSize,
        minReadsPerCytosine=4, cores=64, parallel=TRUE)
    }
  }, error=function(e) { message("  Error: ", e$message); NULL })

  if (is.null(dmrs) || length(dmrs) == 0) {
    return(data.frame(tp=0, fp=0, fn=length(spike_regions),
                      sensitivity=0, precision=NA, f1=NA))
  }

  # Evaluate against ground truth
  hits_real  <- findOverlaps(dmrs, spike_regions)
  hits_false <- findOverlaps(dmrs,
    gaps(spike_regions, end=248956422))

  tp <- length(unique(subjectHits(hits_real)))   # spike regions recovered
  fp <- length(unique(queryHits(hits_false)))     # DMRs outside spike regions
  fn <- length(spike_regions) - tp               # spike regions missed

  sensitivity <- if ((tp+fn) > 0) round(tp/(tp+fn), 3) else NA
  precision   <- if ((tp+fp) > 0) round(tp/(tp+fp), 3) else NA
  f1 <- if (!is.na(sensitivity) && !is.na(precision) &&
            (sensitivity+precision) > 0) {
    round(2*(sensitivity*precision)/(sensitivity+precision), 3)
  } else NA

  data.frame(tp=tp, fp=fp, fn=fn,
             sensitivity=sensitivity,
             precision=precision, f1=f1)
}

results <- list()

for (dm in delta_meths) {
  message("\n=== delta_meth=", dm, " ===")
  treat_file <- file.path(OUT_DIR,
    paste0("aso_vpa_spiked_dm", gsub("\\.", "p", dm), ".rds"))
  aso_vpa_spiked <- readRDS(treat_file)

  for (method in methods) {
    ws_vec <- if (method == "neighbourhood") 100 else window_sizes
    ker_vec <- if (method == "noise_filter") kernels_nf else "NA"

    for (ws in ws_vec) {
      for (ker in ker_vec) {
        for (strict in c(FALSE, TRUE)) {
          mode_lab <- if (strict) "strict" else "loose"
          message("  ", method, " ws=", ws, " ", mode_lab)

          res <- run_and_evaluate(aso_vpa_spiked, aso_ctrl,
                                  method, ws, ker, strict)
          results[[paste(dm, method, ws, ker, mode_lab, sep="_")]] <-
            cbind(data.frame(
              delta_meth=dm, method=method,
              window_size=ws, kernel=ker, mode=mode_lab),
              res)
        }
      }
    }
  }
}

summary_df <- do.call(rbind, results)
write.csv(summary_df,
          file.path(OUT_DIR, "ground_truth_results.csv"),
          row.names=FALSE)

message("\n=== GROUND TRUTH RESULTS ===")
print(summary_df, row.names=FALSE)
message("\nSaved to: ", OUT_DIR, "/ground_truth_results.csv")
