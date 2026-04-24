.libPaths("~/R/library")
library(DMRcaller)
library(ggplot2)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/dmr_benchmark"
CHROM   <- "chr1"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REGION_END <- 10000000

load_group_chr1 <- function(samples) {
  paths <- file.path(COV_DIR,
    paste0(samples, "_", CHROM, ".CpG_report.txt.gz"))
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

chr1_region <- GRanges("chr1", IRanges(1, REGION_END))

make_tiles <- function(binsize) {
  GRanges("chr1", IRanges(
    start = seq(1, REGION_END, by = binsize),
    width = binsize))
}

run_dmrs <- function(treat, ctrl, method, ws,
                     kernel = "triangular",
                     strict = FALSE) {
  regions <- if (method == "bins") make_tiles(ws) else chr1_region

  if (!strict) {
    pval    <- 0.05
    minCpG  <- 1
    minDiff <- 0.1
    minSize <- 1
    test    <- "fisher"
  } else {
    # Strict settings based on DMRcaller CG vignette defaults
    pval    <- 0.01
    minCpG  <- 4
    minDiff <- if (method == "noise_filter") 0.4 else 0.2
    minSize <- 50
    test    <- if (method == "noise_filter") "score" else "fisher"
  }

  tryCatch({
    if (method == "noise_filter") {
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
    } else if (method == "bins") {
      computeDMRs(treat, ctrl,
        regions                 = regions,
        context                 = "CG",
        method                  = "bins",
        binSize                 = ws,
        test                    = test,
        pValueThreshold         = pval,
        minCytosinesCount       = minCpG,
        minProportionDifference = minDiff,
        minGap                  = 0,
        minSize                 = minSize,
        minReadsPerCytosine     = 4,
        cores                   = 20)
    } else {
      computeDMRs(treat, ctrl,
        regions                 = regions,
        context                 = "CG",
        method                  = "neighbourhood",
        test                    = test,
        pValueThreshold         = pval,
        minCytosinesCount       = minCpG,
        minProportionDifference = minDiff,
        minGap                  = 0,
        minSize                 = minSize,
        minReadsPerCytosine     = 4,
        cores                   = 20)
    }
  }, error = function(e) { message("  Error: ", e$message); NULL })
}

results      <- list()
window_sizes <- c(100, 200, 300, 500)
methods      <- c("bins", "neighbourhood", "noise_filter")
kernels_nf   <- c("triangular", "uniform", "epanechnikov")
modes        <- c(FALSE, TRUE)

for (method in methods) {
  for (ws in window_sizes) {
    if (method == "neighbourhood" && ws > 100) next
    if (method == "noise_filter"  && ws > 100) next  # only 100bp for noise_filter
    ker_vec <- if (method == "noise_filter") kernels_nf else "NA"
    for (strict in modes) {
      for (ker in ker_vec) {
        mode_lab <- if (strict) "strict" else "loose"
        message("\n--- ", mode_lab, " | ", method,
                " | ws=", ws, " | kernel=", ker, " ---")

        message("  Real data...")
        dmrs_real <- run_dmrs(aso_vpa, aso_ctrl, method, ws,
                              kernel = ker, strict = strict)
        n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
        message("  Real DMRs: ", n_real)

        message("  Scrambled...")
        dmrs_scr <- run_dmrs(aso_vpa_scr, aso_ctrl_scr, method, ws,
                             kernel = ker, strict = strict)
        n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
        message("  Scrambled DMRs: ", n_scr)

        ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) {
          round(n_real / n_scr, 2)
        } else if (!is.na(n_real) && !is.na(n_scr) && n_scr == 0) {
          Inf
        } else NA
        message("  Signal/Noise: ", ratio)

        results[[paste(method, ws, mode_lab, ker, sep="_")]] <-
          data.frame(
            method      = method,
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
}

summary_df <- do.call(rbind, results)
write.csv(summary_df,
          file.path(OUT_DIR, "parameter_benchmark_chr1_final.csv"),
          row.names = FALSE)

message("\n=== BENCHMARK RESULTS ===")
message(sprintf("%-8s %-15s %-6s %-14s %-12s %-12s %-10s",
                "Mode","Method","ws","Kernel","Real DMRs","Scr DMRs","S/N"))
message(paste(rep("-", 80), collapse=""))
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i,]
  message(sprintf("%-8s %-15s %-6s %-14s %-12s %-12s %-10s",
                  r$mode, r$method, r$window_size, r$kernel,
                  r$n_real, r$n_scrambled, r$ratio))
}

# Counts plot
summary_df$window_size <- factor(summary_df$window_size)
df_long <- reshape(summary_df,
                   varying   = c("n_real","n_scrambled"),
                   v.names   = "n_dmrs",
                   timevar   = "type",
                   times     = c("Real","Scrambled"),
                   direction = "long")

p1 <- ggplot(df_long, aes(x=window_size, y=n_dmrs,
                           colour=type, group=interaction(type,kernel))) +
  geom_line() + geom_point(size=2) +
  facet_grid(mode ~ method) +
  scale_colour_manual(values=c(Real="#02C39A", Scrambled="#F59E0B")) +
  labs(title = "DMR counts: real vs scrambled — chr1 first 10Mb",
       x = "Window/bin size (bp)", y = "Number of DMRs", colour = "") +
  theme_bw(base_size=11)
ggsave(file.path(OUT_DIR, "benchmark_counts_final.pdf"), p1, width=12, height=6)

# S/N ratio plot
summary_df$ratio_num <- as.numeric(ifelse(is.infinite(summary_df$ratio),
                                          NA, summary_df$ratio))
summary_df$method_kernel <- ifelse(summary_df$kernel == "NA",
                                   summary_df$method,
                                   paste0(summary_df$method,"_",summary_df$kernel))

p2 <- ggplot(summary_df, aes(x=window_size, y=ratio_num,
                              colour=method_kernel,
                              group=method_kernel)) +
  geom_line() + geom_point(size=3) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  facet_wrap(~mode) +
  labs(title = "Signal-to-noise ratio (real/scrambled) — chr1 first 10Mb",
       x = "Window/bin size (bp)", y = "Ratio (real/scrambled)",
       colour = "Method/Kernel") +
  theme_bw(base_size=11)
ggsave(file.path(OUT_DIR, "benchmark_ratio_final.pdf"), p2, width=12, height=5)

message("\nDone. Outputs in: ", OUT_DIR)
