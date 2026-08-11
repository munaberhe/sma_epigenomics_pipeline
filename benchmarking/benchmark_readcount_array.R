.libPaths("~/R/library")
library(DMRcaller)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 248956422
N_PERMS    <- 20

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

TASK_ID      <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
window_sizes <- c(100, 200, 300, 500, 1000, 2000)
ws           <- window_sizes[TASK_ID]
message("Task ", TASK_ID, ": window_size=", ws)

scramble_data <- function(dat1, dat2, seed=42) {
  set.seed(seed)
  idx1 <- sample(seq_len(length(dat1)), replace=FALSE)
  idx2 <- sample(seq_len(length(dat2)), replace=FALSE)
  dat1_scr <- dat1; dat2_scr <- dat2
  mcols(dat1_scr)$readsM <- mcols(dat1)$readsM[idx1]
  mcols(dat1_scr)$readsN <- mcols(dat1)$readsN[idx1]
  mcols(dat2_scr)$readsM <- mcols(dat2)$readsM[idx2]
  mcols(dat2_scr)$readsN <- mcols(dat2)$readsN[idx2]
  list(dat1=dat1_scr, dat2=dat2_scr)
}

message("Loading chr1 data...")
aso_vpa <- readBismarkPool(c(
  file.path(COV_DIR, "ASO_VPA_1_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_2_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_3_chr1.CpG_report.txt.gz")))
aso_ctrl <- readBismarkPool(c(
  file.path(COV_DIR, "ASO_CTRL_1_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_CTRL_2_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_CTRL_3_chr1.CpG_report.txt.gz")))
message("Data loaded.")

regions <- GRanges(CHROM, IRanges(start=seq(1, REGION_END, by=ws), width=ws))

message("Real DMR run...")
dmrs_real <- computeDMRs(aso_vpa, aso_ctrl,
  regions=regions, context="CG", method="bins", binSize=ws,
  test="fisher", pValueThreshold=0.01,
  minCytosinesCount=4, minProportionDifference=0.2,
  minGap=0, minSize=50, minReadsPerCytosine=4,
  cores=32, parallel=TRUE)
n_real <- length(dmrs_real)
message("  Real DMRs: ", n_real)

message("20 permutations...")
scr_counts <- integer(N_PERMS)
for (s in 1:N_PERMS) {
  scr <- scramble_data(aso_vpa, aso_ctrl, seed=s)
  d_s <- computeDMRs(scr$dat1, scr$dat2,
    regions=regions, context="CG", method="bins", binSize=ws,
    test="fisher", pValueThreshold=0.01,
    minCytosinesCount=4, minProportionDifference=0.2,
    minGap=0, minSize=50, minReadsPerCytosine=4,
    cores=32, parallel=TRUE)
  scr_counts[s] <- length(d_s)
  message("  seed ", s, "/", N_PERMS, ": ", scr_counts[s])
}

mean_scr <- mean(scr_counts)
sd_scr   <- sd(scr_counts)
ratio    <- ifelse(mean_scr==0, Inf, round(n_real/mean_scr, 3))
message("mean +/- sd: ", round(mean_scr,1), " +/- ", round(sd_scr,1))
message("S/N: ", ratio)

result <- data.frame(
  window_size=ws, n_real=n_real,
  mean_scrambled=mean_scr, sd_scrambled=sd_scr,
  ratio=ratio,
  ratio_lower=round(n_real/(mean_scr+sd_scr), 3),
  ratio_upper=round(n_real/pmax(mean_scr-sd_scr, 0.001), 3),
  n_seeds=N_PERMS,
  method="bins", scramble_method="archie_scramble",
  mode="strict", kernel="NA"
)

out_csv <- file.path(OUT_DIR,
  sprintf("parameter_benchmark_readcount_ws%d.csv", ws))
write.csv(result, out_csv, row.names=FALSE)
message("Saved: ", out_csv)
