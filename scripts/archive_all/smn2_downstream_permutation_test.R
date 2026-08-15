.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

BY_CHR_UNMASK <- "results/alignments/bs/by_chr"
WINDOW <- list(chr="chr5", start=70088223, end=70088522)
N_PERM <- 2000
MIN_READS <- 4

CONDITIONS <- c("ASO_CTRL", "Scramble_CTRL", "ASO_VPA", "Scramble_VPA")

read_unmasked_cpg <- function(condition, chr) {
  files <- file.path(BY_CHR_UNMASK,
                     sprintf("%s_%d_%s.CpG_report.txt.gz", condition, 1:3, chr))
  files <- files[file.exists(files)]
  if (length(files) == 0) return(NULL)
  grs <- lapply(files, function(f) {
    d <- read.table(gzfile(f), header=FALSE, sep="\t",
                    col.names=c("chr","pos","strand","countM","countU","context","tri"),
                    colClasses=c("character","integer","character","integer",
                                 "integer","character","character"))
    d <- d[d$context=="CG", ]
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos,d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

get_proportion <- function(pooled, region) {
  hits <- subsetByOverlaps(pooled, region)
  hits <- hits[hits$readsN >= MIN_READS]
  if (length(hits) == 0) return(list(prop=NA, n=0, readsM=integer(0), readsN=integer(0)))
  list(prop = sum(hits$readsM) / sum(hits$readsN),
       n = length(hits), readsM = hits$readsM, readsN = hits$readsN)
}

message("Loading chr5 methylation for all four conditions...")
region <- GRanges(WINDOW$chr, IRanges(WINDOW$start, WINDOW$end))

obs <- lapply(CONDITIONS, function(cond) {
  message("  ", cond)
  pooled <- read_unmasked_cpg(cond, WINDOW$chr)
  get_proportion(pooled, region)
})
names(obs) <- CONDITIONS

aso_effect <- obs$ASO_CTRL$prop - obs$Scramble_CTRL$prop
vpa_effect <- obs$Scramble_VPA$prop - obs$Scramble_CTRL$prop
predicted  <- obs$Scramble_CTRL$prop + aso_effect + vpa_effect
actual     <- obs$ASO_VPA$prop
deviation  <- actual - predicted

cat("\n=== Observed values ===\n")
cat(sprintf("ASO_CTRL:      %.4f (n=%d CpGs)\n", obs$ASO_CTRL$prop, obs$ASO_CTRL$n))
cat(sprintf("Scramble_CTRL: %.4f (n=%d CpGs)\n", obs$Scramble_CTRL$prop, obs$Scramble_CTRL$n))
cat(sprintf("ASO_VPA:       %.4f (n=%d CpGs)\n", obs$ASO_VPA$prop, obs$ASO_VPA$n))
cat(sprintf("Scramble_VPA:  %.4f (n=%d CpGs)\n", obs$Scramble_VPA$prop, obs$Scramble_VPA$n))
cat(sprintf("\nObserved deviation from additive prediction: %+.4f\n", deviation))

# ---- permutation null ----
message("\nRunning ", N_PERM, " permutations...")
all_M <- c(obs$ASO_CTRL$readsM, obs$Scramble_CTRL$readsM,
           obs$ASO_VPA$readsM, obs$Scramble_VPA$readsM)
all_N <- c(obs$ASO_CTRL$readsN, obs$Scramble_CTRL$readsN,
           obs$ASO_VPA$readsN, obs$Scramble_VPA$readsN)
group_sizes <- c(obs$ASO_CTRL$n, obs$Scramble_CTRL$n, obs$ASO_VPA$n, obs$Scramble_VPA$n)
n_total <- length(all_M)
cuts <- cumsum(group_sizes)

set.seed(42)
perm_devs <- numeric(N_PERM)
for (p in seq_len(N_PERM)) {
  idx <- sample(n_total)
  g1 <- idx[1:cuts[1]]
  g2 <- idx[(cuts[1]+1):cuts[2]]
  g3 <- idx[(cuts[2]+1):cuts[3]]
  g4 <- idx[(cuts[3]+1):cuts[4]]
  p_aso_ctrl <- sum(all_M[g1])/sum(all_N[g1])
  p_scr_ctrl <- sum(all_M[g2])/sum(all_N[g2])
  p_aso_vpa  <- sum(all_M[g3])/sum(all_N[g3])
  p_scr_vpa  <- sum(all_M[g4])/sum(all_N[g4])
  perm_aso <- p_aso_ctrl - p_scr_ctrl
  perm_vpa <- p_scr_vpa - p_scr_ctrl
  perm_pred <- p_scr_ctrl + perm_aso + perm_vpa
  perm_devs[p] <- p_aso_vpa - perm_pred
}

perm_pval <- mean(abs(perm_devs) >= abs(deviation))
cat(sprintf("\n=== Permutation test result (n=%d permutations) ===\n", N_PERM))
cat(sprintf("Permutation p-value: %.4f\n", perm_pval))
cat(sprintf("Permuted deviation distribution: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(perm_devs), sd(perm_devs), min(perm_devs), max(perm_devs)))

if (perm_pval < 0.05) {
  cat("\nRESULT: deviation is larger than expected by chance (p < 0.05)\n")
} else {
  cat("\nRESULT: deviation is consistent with chance given the small sample size at this window (p >= 0.05)\n")
}
