# parameter_benchmark_archie_turbo.R
# PLAN C: archie (read-count permutation) scramble, 32-worker serial.
source("scripts/parameter_benchmark_helpers_turbo.R")
set.seed(42)

dat <- load_chr1_data()
aso_vpa  <- dat$aso_vpa
aso_ctrl <- dat$aso_ctrl

scramble_archie <- function(gr) {
  n <- length(gr)
  perm <- sample.int(n)
  m_ord <- mcols(gr)
  if ("readsM" %in% colnames(m_ord)) {
    m_ord$readsM <- m_ord$readsM[perm]
    m_ord$readsN <- m_ord$readsN[perm]
  }
  mcols(gr) <- m_ord
  gr
}

aso_vpa_scr  <- scramble_archie(aso_vpa)
aso_ctrl_scr <- scramble_archie(aso_ctrl)

run_benchmark_loop(aso_vpa, aso_ctrl, aso_vpa_scr, aso_ctrl_scr, null_name = "archie_scramble")
