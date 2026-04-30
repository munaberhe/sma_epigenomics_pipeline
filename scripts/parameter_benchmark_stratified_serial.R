.libPaths("~/R/library")
# parameter_benchmark_stratified_serial.R — PLAN B (stratified scramble)
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

.this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) > 0) dirname(normalizePath(f)) else getwd()
  }
)
source(file.path(.this_dir, "parameter_benchmark_helpers_serial.R"))

NULL_NAME <- "stratified_scramble"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dat <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl

scramble_stratified <- function(dat1, dat2, seed=42) {
  set.seed(seed)
  coverage <- mcols(dat1)$readsN
  strata <- cut(coverage, breaks=c(0,5,10,20,50,Inf),
                labels=c("1-5","6-10","11-20","21-50","50+"),
                include.lowest=TRUE)
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
  message("  Strata sizes: ", paste(table(strata), collapse=", "))
  message("  Original  ASO_VPA mean methylation:  ",
          round(sum(mcols(dat1)$readsM)/sum(mcols(dat1)$readsN)*100,1), "%")
  message("  Scrambled ASO_VPA mean methylation:  ",
          round(sum(mcols(dat1_scr)$readsM)/sum(mcols(dat1_scr)$readsN)*100,1), "%")
  message("  Original  ASO_CTRL mean methylation: ",
          round(sum(mcols(dat2)$readsM)/sum(mcols(dat2)$readsN)*100,1), "%")
  message("  Scrambled ASO_CTRL mean methylation: ",
          round(sum(mcols(dat2_scr)$readsM)/sum(mcols(dat2_scr)$readsN)*100,1), "%")
  list(dat1=dat1_scr, dat2=dat2_scr)
}

message("Creating stratified scrambled datasets...")
scrambled <- scramble_stratified(aso_vpa, aso_ctrl, seed=42)
aso_vpa_scr  <- scrambled$dat1
aso_ctrl_scr <- scrambled$dat2
message("  Done")

run_benchmark_loop(aso_vpa, aso_ctrl, aso_vpa_scr, aso_ctrl_scr, NULL_NAME)
