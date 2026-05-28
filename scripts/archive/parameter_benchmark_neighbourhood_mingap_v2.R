.libPaths("~/R/library")
# parameter_benchmark_neighbourhood_mingap_v2.R
# Adds archie + stratified null models to existing label_swap mingap results
# Uses same parameters as original run (minDiff=0.2 strict, minDiff=0.1 loose)
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(BiocParallel)

COV_DIR  <- "results/alignments/bs/by_chr"
OUT_DIR  <- "results/dmr_benchmark"
CHROM    <- "chr1"
RDS_PATH <- "results/dmr_benchmark/chr1_data.rds"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REGION_END    <- 248956422
MINGAP_VALUES <- c(100, 200, 500, 1000, 2000)
BPPARAM       <- SnowParam(workers = 16, type = "SOCK")

message("Loading data from cache...")
dat      <- readRDS(RDS_PATH)
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl
message("  Loaded.")

scramble_stratified <- function(dat1, dat2, seed=42) {
  set.seed(seed)
  coverage <- mcols(dat1)$readsN
  strata <- cut(coverage, breaks=c(0,5,10,20,50,Inf),
                labels=c("1-5","6-10","11-20","21-50","50+"),
                include.lowest=TRUE)
  idx <- seq_along(coverage)
  for (s in levels(strata)) {
    si <- which(strata == s)
    if (length(si) > 1) idx[si] <- sample(si)
  }
  dat1_scr <- dat1; dat2_scr <- dat2
  mcols(dat1_scr)$readsM <- mcols(dat1)$readsM[idx]
  mcols(dat1_scr)$readsN <- mcols(dat1)$readsN[idx]
  mcols(dat2_scr)$readsM <- mcols(dat2)$readsM[idx]
  mcols(dat2_scr)$readsN <- mcols(dat2)$readsN[idx]
  list(dat1=dat1_scr, dat2=dat2_scr)
}

scramble_counts <- function(dat, seed=42) {
  set.seed(seed)
  dat_scr <- dat
  idx <- sample(length(dat))
  mcols(dat_scr)$readsM <- mcols(dat)$readsM[idx]
  mcols(dat_scr)$readsN <- mcols(dat)$readsN[idx]
  dat_scr
}

message("Creating scrambled datasets...")
strat          <- scramble_stratified(aso_vpa, aso_ctrl, seed=42)
aso_vpa_strat  <- strat$dat1
aso_ctrl_strat <- strat$dat2
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

ckpt    <- file.path(OUT_DIR, "checkpoint_neighbourhood_mingap_v2.rds")
results <- if (file.exists(ckpt)) { message("Resuming from checkpoint..."); readRDS(ckpt) } else list()

null_models <- list(
  stratified = list(vpa=aso_vpa_strat,  ctrl=aso_ctrl_strat),
  archie     = list(vpa=aso_vpa_archie, ctrl=aso_ctrl_archie)
)

for (nm in names(null_models)) {
  scr_vpa  <- null_models[[nm]]$vpa
  scr_ctrl <- null_models[[nm]]$ctrl
  for (mg in MINGAP_VALUES) {
    for (mode in c("loose","strict")) {
      key <- paste(nm, mg, mode, sep="_")
      if (!is.null(results[[key]])) { message("skip: ", key); next }
      message("\n--- ", nm, " | minGap=", mg, " | ", mode, " ---")
      message("  Real...")
      dmrs_real <- run_neighbourhood(aso_vpa, aso_ctrl, mg, mode)
      n_real    <- if (!is.null(dmrs_real)) length(dmrs_real) else NA
      message("  Real DMRs: ", n_real)
      message("  Scrambled...")
      dmrs_scr <- run_neighbourhood(scr_vpa, scr_ctrl, mg, mode)
      n_scr    <- if (!is.null(dmrs_scr)) length(dmrs_scr) else NA
      message("  Scrambled DMRs: ", n_scr)
      ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) round(n_real/n_scr,2) else
               if (!is.na(n_scr) && n_scr == 0) Inf else NA
      message("  S/N: ", ratio)
      results[[key]] <- data.frame(
        method="neighbourhood_mingap", window_size=mg, mode=mode, kernel="NA",
        n_real=n_real, n_scrambled=n_scr, ratio=ratio,
        scramble_method=nm, stringsAsFactors=FALSE)
      saveRDS(results, ckpt)
    }
  }
}

# Combine with existing label_swap CSV
existing <- read.csv(file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap.csv"))
if (!"scramble_method" %in% colnames(existing)) existing$scramble_method <- "label_swap"
new_df   <- do.call(rbind, results)
combined <- rbind(existing, new_df)

write.csv(combined,
          file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap_combined.csv"),
          row.names=FALSE)
message("\nDone. Saved: parameter_benchmark_neighbourhood_mingap_combined.csv")
