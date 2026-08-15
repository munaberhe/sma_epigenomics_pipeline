.libPaths("~/R/library")
# nb_mingap_archie_only.R
# Archie scramble only for neighbourhood minGap benchmark
# Writes to separate checkpoint so it doesn't conflict with v2 job
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(BiocParallel)

COV_DIR  <- "results/alignments/bs/by_chr"
OUT_DIR  <- "results/dmr_benchmark"
CHROM    <- "chr1"
RDS_PATH <- "results/dmr_benchmark/chr1_data.rds"

REGION_END    <- 248956422
MINGAP_VALUES <- c(100, 200, 500, 1000, 2000)
MINSIZE_VALUES <- c(1)
BPPARAM       <- SnowParam(workers = 16, type = "SOCK")

message("Loading data from cache...")
dat      <- readRDS(RDS_PATH)
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl
message("  Loaded.")

scramble_counts <- function(dat, seed=42) {
  set.seed(seed)
  dat_scr <- dat
  idx <- sample(length(dat))
  mcols(dat_scr)$readsM <- mcols(dat)$readsM[idx]
  mcols(dat_scr)$readsN <- mcols(dat)$readsN[idx]
  dat_scr
}

message("Creating archie scrambled datasets...")
aso_vpa_archie  <- scramble_counts(aso_vpa,  seed=42)
aso_ctrl_archie <- scramble_counts(aso_ctrl, seed=123)
message("  Done.")

get_thresholds <- function(mode) {
  if (mode == "loose")
    list(pval=0.05, minCpG=1, minDiff=0.1, minSize=1,  test="fisher")
  else
    list(pval=0.01, minCpG=4, minDiff=0.2, minSize=50, test="fisher")
}

run_neighbourhood <- function(treat, ctrl, mingap, mode) {
  th     <- get_thresholds(mode)
  region <- GRanges(CHROM, IRanges(1, REGION_END))
  tryCatch({
    computeDMRs(treat, ctrl, regions=region, context="CG",
      method="neighbourhood", test=th$test,
      pValueThreshold=th$pval, minCytosinesCount=th$minCpG,
      minProportionDifference=th$minDiff,
      minGap=mingap, minSize=th$minSize,
      minReadsPerCytosine=4, BPPARAM=BPPARAM, parallel=TRUE)
  }, error=function(e) { message("  Error: ", e$message); NULL })
}

ckpt    <- file.path(OUT_DIR, "checkpoint_neighbourhood_mingap_archie.rds")
results <- if (file.exists(ckpt)) { message("Resuming..."); readRDS(ckpt) } else list()

for (mg in MINGAP_VALUES) {
  for (mode in c("loose","strict")) {
    key <- paste("archie", mg, mode, sep="_")
    if (!is.null(results[[key]])) { message("skip: ", key); next }
    message("\n--- archie | minGap=", mg, " | ", mode, " ---")
    message("  Real...")
    dmrs_real <- run_neighbourhood(aso_vpa, aso_ctrl, mg, mode)
    n_real    <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
    message("  Real DMRs: ", n_real)
    message("  Scrambled...")
    dmrs_scr <- run_neighbourhood(aso_vpa_archie, aso_ctrl_archie, mg, mode)
    n_scr    <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
    message("  Scrambled DMRs: ", n_scr)
    ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) round(n_real/n_scr,2) else
             if (!is.na(n_scr) && n_scr == 0) Inf else NA
    message("  S/N: ", ratio)
    results[[key]] <- data.frame(
      method="neighbourhood_mingap", window_size=mg, mode=mode, kernel="NA",
      n_real=n_real, n_scrambled=n_scr, ratio=ratio,
      scramble_method="archie", stringsAsFactors=FALSE)
    saveRDS(results, ckpt)
  }
}

summary_df <- do.call(rbind, results)
write.csv(summary_df,
          file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap_archie.csv"),
          row.names=FALSE)
message("\nDone. Saved: parameter_benchmark_neighbourhood_mingap_archie.csv")
