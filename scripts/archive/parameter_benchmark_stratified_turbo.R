# parameter_benchmark_stratified_turbo.R
# PLAN C: stratified (coverage-bin) scramble, 32-worker serial.
source("scripts/parameter_benchmark_helpers_turbo.R")
set.seed(42)

dat <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl

scramble_stratified <- function(gr, n_bins = 10) {
  m_ord <- mcols(gr)
  if (!("readsN" %in% colnames(m_ord))) return(gr)
  cov <- m_ord$readsN
  bins <- cut(cov,
              breaks = unique(quantile(cov, probs = seq(0, 1, length.out = n_bins + 1),
                                        na.rm = TRUE)),
              include.lowest = TRUE, labels = FALSE)
  perm <- seq_along(cov)
  for (b in unique(bins[!is.na(bins)])) {
    idx <- which(bins == b)
    if (length(idx) > 1) perm[idx] <- sample(idx)
  }
  m_ord$readsM <- m_ord$readsM[perm]
  m_ord$readsN <- m_ord$readsN[perm]
  mcols(gr) <- m_ord
  gr
}

aso_vpa_scr  <- scramble_stratified(aso_vpa)
aso_ctrl_scr <- scramble_stratified(aso_ctrl)

run_benchmark_loop(aso_vpa, aso_ctrl, aso_vpa_scr, aso_ctrl_scr, null_name = "stratified_scramble")
