.libPaths("~/R/library")
# parameter_benchmark_neighbourhood_mingap_final.R
# Neighbourhood minGap benchmark using Radu's exact parameters
# minGap varied at 100, 200, 300, 500, 1000, 2000bp on whole chr1
# minSize tested at 100bp and 200bp as per Radu's suggestion
# Three null models: label swap, Archie scramble, stratified scramble
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

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

# --- Null model 1: Label swap ---
aso_vpa_ls  <- aso_ctrl
aso_ctrl_ls <- aso_vpa

# --- Null model 2: Archie scramble ---
scramble_archie <- function(dat, seed=42) {
  set.seed(seed)
  dat_scr <- dat
  idx <- sample(length(dat))
  mcols(dat_scr)$readsM <- mcols(dat)$readsM[idx]
  mcols(dat_scr)$readsN <- mcols(dat)$readsN[idx]
  dat_scr
}
aso_vpa_ar  <- scramble_archie(aso_vpa,  seed=42)
aso_ctrl_ar <- scramble_archie(aso_ctrl, seed=123)

# --- Null model 3: Stratified scramble ---
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
    if (length(stratum_idx) > 1) idx[stratum_idx] <- sample(stratum_idx)
  }
  dat1_scr <- dat1; dat2_scr <- dat2
  mcols(dat1_scr)$readsM <- mcols(dat1)$readsM[idx]
  mcols(dat1_scr)$readsN <- mcols(dat1)$readsN[idx]
  mcols(dat2_scr)$readsM <- mcols(dat2)$readsM[idx]
  mcols(dat2_scr)$readsN <- mcols(dat2)$readsN[idx]
  list(dat1=dat1_scr, dat2=dat2_scr)
}
message("Creating scrambled datasets...")
strat <- scramble_stratified(aso_vpa, aso_ctrl, seed=42)
aso_vpa_st  <- strat$dat1
aso_ctrl_st <- strat$dat2
message("  Done")

run_neighbourhood <- function(treat, ctrl, mingap, minsize) {
  tryCatch({
    computeDMRs(treat, ctrl,
      regions                 = chr1_region,
      context                 = "CG",
      method                  = "neighbourhood",
      test                    = "score",
      pValueThreshold         = 0.01,
      minCytosinesCount       = 4,
      minProportionDifference = 0.4,
      minGap                  = mingap,
      minSize                 = minsize,
      minReadsPerCytosine     = 4,
      cores                   = 32,
      parallel                = TRUE)
  }, error = function(e) { message("  Error: ", e$message); NULL })
}

mingap_sizes <- c(100, 200, 300, 500, 1000, 2000)
minsize_vals <- c(100, 200)
null_models  <- list(
  label_swap = list(treat=aso_vpa_ls,  ctrl=aso_ctrl_ls),
  archie     = list(treat=aso_vpa_ar,  ctrl=aso_ctrl_ar),
  stratified = list(treat=aso_vpa_st,  ctrl=aso_ctrl_st)
)
results <- list()

for (nm in names(null_models)) {
  treat <- null_models[[nm]]$treat
  ctrl  <- null_models[[nm]]$ctrl
  for (mg in mingap_sizes) {
    for (ms in minsize_vals) {
      message("\n--- ", nm, " | minGap=", mg, " | minSize=", ms, " ---")

      message("  Real data...")
      dmrs_real <- run_neighbourhood(aso_vpa, aso_ctrl, mg, ms)
      n_real <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
      message("  Real DMRs: ", n_real)

      message("  Scrambled...")
      dmrs_scr <- run_neighbourhood(treat, ctrl, mg, ms)
      n_scr <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
      message("  Scrambled DMRs: ", n_scr)

      ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) {
        round(n_real / n_scr, 2)
      } else if (!is.na(n_real) && !is.na(n_scr) && n_scr == 0) {
        Inf
      } else NA
      message("  Signal/Noise: ", ratio)

      results[[paste("neighbourhood", nm, mg, ms, sep="_")]] <- data.frame(
        method          = "neighbourhood_mingap",
        mingap          = mg,
        minsize         = ms,
        scramble_method = nm,
        n_real          = n_real,
        n_scrambled     = n_scr,
        ratio           = ratio,
        stringsAsFactors = FALSE
      )
    }
  }
}

summary_df <- do.call(rbind, results)
write.csv(summary_df,
          file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap_final.csv"),
          row.names = FALSE)

message("\n=== NEIGHBOURHOOD minGap RESULTS ===")
message(sprintf("%-12s %-6s %-8s %-12s %-12s %-10s",
                "Null model", "minGap", "minSize", "Real DMRs", "Scr DMRs", "S/N"))
message(paste(rep("-", 65), collapse=""))
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i,]
  message(sprintf("%-12s %-6s %-8s %-12s %-12s %-10s",
                  r$scramble_method, r$mingap, r$minsize,
                  r$n_real, r$n_scrambled, r$ratio))
}
message("\nDone. Saved to: ", OUT_DIR,
        "/parameter_benchmark_neighbourhood_mingap_final.csv")
