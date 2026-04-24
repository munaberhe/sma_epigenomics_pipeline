
.libPaths("~/R/library")
# benchmark_parameters.R
# Parameter optimisation for DMR calling
# Tests noise_filter, bins, neighbourhood across window sizes 100, 200, 300, 500bp
# Real data vs scrambled data on chr1 — following Zabet et al. DMRcaller paper
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

COV_DIR  <- "results/alignments/bs/by_chr"
OUT_DIR  <- "results/dmr_benchmark"
CHROM    <- "chr1"
CONTEXT  <- "CG"
MIN_COV  <- 4
SET.SEED <- 42

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Use ASO_VPA vs ASO_CTRL as the key contrast
TREATMENT_SAMPLES <- c("ASO_VPA_1", "ASO_VPA_2", "ASO_VPA_3")
CONTROL_SAMPLES   <- c("ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3")

# Load chr1 data
message("Loading chr1 data...")
load_group <- function(samples) {
  paths <- file.path(COV_DIR,
    paste0(samples, "_", CHROM, ".CpG_report.txt.gz"))
  paths <- paths[file.exists(paths)]
  dat <- readBismarkPool(paths)
  dat[dat$readsN >= MIN_COV]
}

treatment_real <- load_group(TREATMENT_SAMPLES)
control_real   <- load_group(CONTROL_SAMPLES)
message("  Treatment CpGs: ", length(treatment_real))
message("  Control CpGs:   ", length(control_real))

# Scramble function — shuffle readsM values across positions
# Preserves coverage (readsN) but destroys spatial methylation pattern
scramble_methylation <- function(methData, seed = SET.SEED) {
  set.seed(seed)
  scrambled <- methData
  idx <- sample(length(scrambled))
  scrambled$readsM <- methData$readsM[idx]
  scrambled
}

message("Creating scrambled datasets...")
treatment_scr <- scramble_methylation(treatment_real, seed = 42)
control_scr   <- scramble_methylation(control_real,   seed = 123)
message("  Scrambled datasets created")

# Parameter grid
methods      <- c("noise_filter", "bins", "neighbourhood")
window_sizes <- c(100, 200, 300, 500)

# Results storage
results <- data.frame()

for (method in methods) {
  for (ws in window_sizes) {

    # neighbourhood method doesn't use windowSize — skip 200/300/500
    if (method == "neighbourhood" && ws != 100) next

    message("\n--- Method: ", method, " | Window/Bin: ", ws, "bp ---")

    run_dmrs <- function(t_dat, c_dat, label) {
      tryCatch({
        if (method == "noise_filter") {
          computeDMRs(t_dat, c_dat,
            regions                 = GRanges(seqnames=Rle(CHROM),
                                              ranges=IRanges(1, max(end(t_dat)))),
            context                 = CONTEXT,
            method                  = "noise_filter",
            windowSize              = ws,
            kernelFunction          = "triangular",
            test                    = "score",
            pValueThreshold         = 0.01,
            minCytosinesCount       = 4,
            minProportionDifference = 0.1,
            minGap                  = 0,
            minSize                 = 50,
            minReadsPerCytosine     = 4,
            cores                   = 4)
        } else if (method == "bins") {
          computeDMRs(t_dat, c_dat,
            regions                 = GRanges(seqnames=Rle(CHROM),
                                              ranges=IRanges(1, max(end(t_dat)))),
            context                 = CONTEXT,
            method                  = "bins",
            binSize                 = ws,
            test                    = "score",
            pValueThreshold         = 0.01,
            minCytosinesCount       = 4,
            minProportionDifference = 0.1,
            minGap                  = 200,
            minSize                 = 50,
            minReadsPerCytosine     = 4,
            cores                   = 4)
        } else {
          computeDMRs(t_dat, c_dat,
            regions                 = GRanges(seqnames=Rle(CHROM),
                                              ranges=IRanges(1, max(end(t_dat)))),
            context                 = CONTEXT,
            method                  = "neighbourhood",
            test                    = "score",
            pValueThreshold         = 0.01,
            minCytosinesCount       = 4,
            minProportionDifference = 0.1,
            minGap                  = 200,
            minSize                 = 1,
            minReadsPerCytosine     = 4,
            cores                   = 4)
        }
      }, error = function(e) {
        message("  Error: ", e$message)
        NULL
      })
    }

    # Real data
    message("  Running on real data...")
    dmrs_real <- run_dmrs(treatment_real, control_real, "real")
    n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
    message("  Real DMRs: ", n_real)

    # Scrambled data
    message("  Running on scrambled data...")
    dmrs_scr <- run_dmrs(treatment_scr, control_scr, "scrambled")
    n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
    message("  Scrambled DMRs: ", n_scr)

    # Signal to noise ratio
    snr <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) {
      round(n_real / n_scr, 2)
    } else if (!is.na(n_real) && n_scr == 0) {
      Inf
    } else NA

    results <- rbind(results, data.frame(
      method      = method,
      window_size = ws,
      real_dmrs   = n_real,
      scr_dmrs    = n_scr,
      signal_noise_ratio = snr,
      stringsAsFactors = FALSE
    ))

    message("  Signal/Noise ratio: ", snr)
  }
}

# Save results
write.csv(results, file.path(OUT_DIR, "benchmark_results.csv"), row.names = FALSE)

# Print summary table
message("\n=== BENCHMARK RESULTS ===")
message(sprintf("%-15s %-12s %-12s %-12s %-10s",
                "Method", "Window(bp)", "Real DMRs", "Scr DMRs", "S/N Ratio"))
message(paste(rep("-", 65), collapse=""))
for (i in seq_len(nrow(results))) {
  r <- results[i,]
  message(sprintf("%-15s %-12s %-12s %-12s %-10s",
                  r$method, r$window_size,
                  r$real_dmrs, r$scr_dmrs,
                  r$signal_noise_ratio))
}

# Plot benchmark results
pdf(file.path(OUT_DIR, "benchmark_plot.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

cols <- c("noise_filter"="#02C39A", "bins"="#F59E0B", "neighbourhood"="#065A82")

# Plot real DMRs
plot(NULL, xlim=c(50,550), ylim=c(0, max(results$real_dmrs, na.rm=TRUE)*1.1),
     xlab="Window/Bin size (bp)", ylab="Number of DMRs",
     main="Real Data — DMR counts by method and window size")
for (m in methods) {
  sub <- results[results$method == m, ]
  lines(sub$window_size, sub$real_dmrs, col=cols[m], lwd=2, pch=16, type="b")
}
legend("topright", legend=methods, col=cols, lwd=2, bty="n")

# Plot signal/noise
plot(NULL, xlim=c(50,550),
     ylim=c(0, max(results$signal_noise_ratio[is.finite(results$signal_noise_ratio)],
                   na.rm=TRUE)*1.1),
     xlab="Window/Bin size (bp)", ylab="Signal/Noise ratio (real/scrambled)",
     main="Signal to Noise Ratio")
for (m in methods) {
  sub <- results[results$method == m, ]
  sub$snr <- ifelse(is.infinite(sub$signal_noise_ratio), NA, sub$signal_noise_ratio)
  lines(sub$window_size, sub$snr, col=cols[m], lwd=2, pch=16, type="b")
}
legend("topright", legend=methods, col=cols, lwd=2, bty="n")
dev.off()

message("\nBenchmark complete. Results in: ", OUT_DIR)
