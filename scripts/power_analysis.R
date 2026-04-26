.libPaths("~/R/library")
# power_analysis.R — fixed binomial model
# Each CpG has readsN total reads, readsM methylated
# Two groups compared: n1 reads in group1, n2 reads in group2

min_detectable <- function(n1, n2=n1, p1=0.64, alpha=0.01, power=0.80) {
  for (d in seq(0.001, 0.60, by=0.001)) {
    p2 <- p1 - d
    if (p2 < 0) break
    # Pooled proportion under H0
    p_pool <- (n1*p1 + n2*p2) / (n1 + n2)
    # SE of difference under H0 (binomial)
    se_h0 <- sqrt(p_pool*(1-p_pool)*(1/n1 + 1/n2))
    # SE under H1
    se_h1 <- sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)
    z_alpha <- qnorm(1 - alpha/2)
    # Power
    achieved <- pnorm((abs(p1-p2) - z_alpha*se_h0) / se_h1)
    if (achieved >= power) return(round(d*100, 2))
  }
  return(NA)
}

cat("\n=== COVERAGE POWER ANALYSIS ===\n")
cat("Binomial model — each CpG has n reads, comparing two groups\n")
cat("alpha=0.01, 80% power, baseline methylation=64%\n\n")

cat("Per-CpG reads (single replicate) vs minimum detectable diff:\n")
for (n in c(4, 5, 8, 10, 12, 15, 20, 30)) {
  mdd <- min_detectable(n, n)
  mdd_val <- ifelse(is.na(mdd), ">60", sprintf("%.1f%%", mdd))
  flag <- ifelse(is.na(mdd) || mdd > 10, " << underpowered", " OK")
  cat(sprintf("  %2dx reads: %s min detectable%s\n", n, mdd_val, flag))
}

cat("\nWith pooled replicates (sum of reads across 3 replicates per group):\n")
for (per_rep in c(4, 10, 12, 14)) {
  n_pooled <- per_rep * 3
  mdd <- min_detectable(n_pooled, n_pooled)
  mdd_val <- ifelse(is.na(mdd), ">60", sprintf("%.1f%%", mdd))
  flag <- ifelse(is.na(mdd) || mdd > 10, " << underpowered", " OK")
  cat(sprintf("  %2dx per replicate (%2dx pooled): %s min detectable%s\n",
    per_rep, n_pooled, mdd_val, flag))
}

cat("\nConclusion:\n")
mdd_4x  <- min_detectable(4, 4)
mdd_10x <- min_detectable(10, 10)
mdd_12x <- min_detectable(12*3, 12*3)
cat(sprintf("  Single replicate 4x:      %.1f%% min detectable\n", mdd_4x))
cat(sprintf("  Single replicate 10x:     %.1f%% min detectable\n", mdd_10x))
cat(sprintf("  Pooled 3 reps at 12x each: %.1f%% min detectable\n", mdd_12x))
