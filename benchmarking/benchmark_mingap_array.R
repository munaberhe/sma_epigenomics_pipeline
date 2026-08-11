.libPaths("~/R/library")
library(DMRcaller)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark"
CHROM      <- "chr1"
REGION_END <- 248956422
N_SEEDS    <- 20

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# Get task parameters from SLURM array index
TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
mingap_values  <- c(100, 200, 300, 500, 1000, 2000)
minsize_values <- c(100, 200)
combos <- expand.grid(mingap=mingap_values, minsize=minsize_values)
mg <- combos$mingap[TASK_ID]
ms <- combos$minsize[TASK_ID]
message("Task ", TASK_ID, ": minGap=", mg, " minSize=", ms)

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

region <- GRanges(CHROM, IRanges(1, REGION_END))

message("Real DMR run...")
dmrs_real <- computeDMRs(aso_vpa, aso_ctrl,
  regions=region, context="CG", method="neighbourhood",
  test="score", pValueThreshold=0.01,
  minCytosinesCount=4, minProportionDifference=0.2,
  minGap=mg, minSize=ms,
  minReadsPerCytosine=4, cores=32, parallel=TRUE)
n_real <- length(dmrs_real)
message("  Real DMRs: ", n_real)

message("20 permutations...")
scr_counts <- integer(N_SEEDS)
for (s in seq_len(N_SEEDS)) {
  scr_s <- scramble_data(aso_vpa, aso_ctrl, seed=s)
  d_s <- computeDMRs(scr_s$dat1, scr_s$dat2,
    regions=region, context="CG", method="neighbourhood",
    test="score", pValueThreshold=0.01,
    minCytosinesCount=4, minProportionDifference=0.2,
    minGap=mg, minSize=ms,
    minReadsPerCytosine=4, cores=32, parallel=TRUE)
  scr_counts[s] <- length(d_s)
  message("  seed ", s, "/", N_SEEDS, ": ", scr_counts[s])
}

mean_scr <- mean(scr_counts)
sd_scr   <- sd(scr_counts)
ratio    <- ifelse(mean_scr==0, Inf, round(n_real/mean_scr, 3))
message("mean +/- sd: ", round(mean_scr,1), " +/- ", round(sd_scr,1))
message("S/N: ", ratio)

result <- data.frame(
  method="neighbourhood", window_size=mg, mode="strict",
  kernel="NA", scramble_method="archie_scramble",
  n_real=n_real, n_scrambled=round(mean_scr),
  ratio=ratio, mean_scrambled=mean_scr,
  sd_scrambled=sd_scr,
  size_real_median=median(width(dmrs_real)),
  size_real_mean=mean(width(dmrs_real)),
  stringsAsFactors=FALSE)

out_csv <- file.path(OUT_DIR,
  sprintf("parameter_benchmark_mingap_mg%d_ms%d.csv", mg, ms))
write.csv(result, out_csv, row.names=FALSE)
message("Saved: ", out_csv)
