.libPaths("~/R/library")
# coverage_power_analysis.R
# Real coverage distribution across all samples and chromosomes
# Calculates what proportion of CpGs are adequately powered to detect
# methylation differences at our threshold
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(ggplot2)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/qc/coverage_power"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

SAMPLES <- c("Scramble_CTRL_1","Scramble_CTRL_2","Scramble_CTRL_3",
             "Scramble_VPA_1", "Scramble_VPA_2", "Scramble_VPA_3",
             "ASO_CTRL_1",     "ASO_CTRL_2",     "ASO_CTRL_3",
             "ASO_VPA_1",      "ASO_VPA_2",      "ASO_VPA_3")

CHROMOSOMES <- paste0("chr", c(1:22, "X"))
BREAKS      <- c(1, 4, 5, 10, 15, 20)

# Minimum detectable difference at given coverage and power
min_detectable_diff <- function(n, p_baseline=0.64, alpha=0.01, power=0.80) {
  for (d in seq(0.01, 0.60, by=0.005)) {
    p2 <- p_baseline - d
    if (p2 < 0) break
    p_pool <- (p_baseline + p2) / 2
    se <- sqrt(2 * p_pool * (1-p_pool) / n)
    z_alpha <- qnorm(1 - alpha/2)
    achieved <- pnorm(abs(p_baseline-p2)/se - z_alpha)
    if (achieved >= power) return(d)
  }
  return(NA)
}

# Accumulate coverage counts across all samples and chromosomes
message("Computing coverage distribution across all samples and chromosomes...")
total_sites <- 0
coverage_counts <- integer(101)  # counts for 0–100x

for (sample in SAMPLES) {
  message("  Sample: ", sample)
  for (chrom in CHROMOSOMES) {
    path <- file.path(COV_DIR,
      paste0(sample, "_", chrom, ".CpG_report.txt.gz"))
    if (!file.exists(path)) next
    tryCatch({
      dat <- readBismark(path)
      dat <- dat[dat$context == "CG"]
      n   <- dat$readsN
      # Bin coverage — cap at 100
      n_capped <- pmin(n, 100)
      for (i in seq_along(n_capped)) {
        coverage_counts[n_capped[i]+1] <-
          coverage_counts[n_capped[i]+1] + 1L
      }
      total_sites <- total_sites + length(n)
    }, error=function(e) message("    Error: ", e$message))
  }
}

message("Total CpG observations: ", total_sites)

# Proportion of sites at each coverage threshold
cov_df <- data.frame(
  coverage     = 0:100,
  count        = coverage_counts,
  proportion   = coverage_counts / total_sites
)

# Cumulative proportion at >= each threshold
cov_summary <- data.frame(
  threshold = BREAKS,
  pct_sites = sapply(BREAKS, function(b)
    sum(coverage_counts[(b+1):101]) / total_sites * 100),
  min_detectable_diff_80pct = sapply(BREAKS, function(b)
    min_detectable_diff(b, power=0.80) * 100),
  min_detectable_diff_90pct = sapply(BREAKS, function(b)
    min_detectable_diff(b, power=0.90) * 100)
)

write.csv(cov_df,      file.path(OUT_DIR, "coverage_distribution_all.csv"),  row.names=FALSE)
write.csv(cov_summary, file.path(OUT_DIR, "coverage_power_summary.csv"), row.names=FALSE)

# Print summary table
message("\n=== COVERAGE POWER SUMMARY (all samples, all chromosomes) ===")
message("Baseline methylation: 64% (controls), alpha=0.01, score test\n")
message(sprintf("%-12s %-18s %-22s %-22s",
  "Min reads", "% CpGs covered", "Min detectable (80%)", "Min detectable (90%)"))
message(paste(rep("-", 76), collapse=""))
for (i in seq_len(nrow(cov_summary))) {
  r <- cov_summary[i,]
  message(sprintf("%-12d %-18s %-22s %-22s",
    r$threshold,
    sprintf("%.1f%%", r$pct_sites),
    sprintf("%.1f%%", r$min_detectable_diff_80pct),
    sprintf("%.1f%%", r$min_detectable_diff_90pct)))
}

message("\n=== KEY FINDINGS FOR RADU ===")
r4  <- cov_summary[cov_summary$threshold==4,  "min_detectable_diff_80pct"]
r10 <- cov_summary[cov_summary$threshold==10, "min_detectable_diff_80pct"]
p4  <- cov_summary[cov_summary$threshold==4,  "pct_sites"]
p10 <- cov_summary[cov_summary$threshold==10, "pct_sites"]
message(sprintf("At minReadsPerCytosine=4:  %.1f%% of CpGs covered, min detectable diff = %.1f%%",
  p4, r4))
message(sprintf("At minReadsPerCytosine=10: %.1f%% of CpGs covered, min detectable diff = %.1f%%",
  p10, r10))
message(sprintf("Our threshold (0.1 = 10%%): %s at 4x coverage",
  ifelse(r4 > 10, "UNDERPOWERED", "adequately powered")))
message(sprintf("Our threshold (0.1 = 10%%): %s at 10x coverage",
  ifelse(r10 > 10, "UNDERPOWERED", "adequately powered")))

# Plot coverage distribution
p1 <- ggplot(cov_df[cov_df$coverage <= 30,], aes(x=coverage, y=proportion*100)) +
  geom_col(fill="#065A82", alpha=0.8) +
  geom_vline(xintercept=4,  colour="#8892A4", linewidth=1, linetype="dashed") +
  geom_vline(xintercept=10, colour="#02C39A", linewidth=1, linetype="dashed") +
  annotate("text", x=4.5,  y=max(cov_df$proportion[1:31]*100)*0.9,
    label="4x filter", colour="#8892A4", hjust=0, size=3.5) +
  annotate("text", x=10.5, y=max(cov_df$proportion[1:31]*100)*0.8,
    label="10x filter", colour="#02C39A", hjust=0, size=3.5) +
  labs(title="CpG coverage distribution — all 12 samples, all chromosomes",
       x="Reads per cytosine", y="% of CpG sites") +
  theme_bw(base_size=12)
ggsave(file.path(OUT_DIR, "coverage_distribution.pdf"), p1, width=10, height=5)

# Plot power curve
power_curve <- data.frame(
  coverage = 1:50,
  min_diff = sapply(1:50, function(n)
    min_detectable_diff(n, power=0.80) * 100)
)
p2 <- ggplot(power_curve, aes(x=coverage, y=min_diff)) +
  geom_line(colour="#065A82", linewidth=1.2) +
  geom_hline(yintercept=10, colour="#F59E0B", linewidth=1, linetype="dotted") +
  geom_vline(xintercept=4,  colour="#8892A4", linewidth=0.8, linetype="dashed") +
  geom_vline(xintercept=10, colour="#02C39A", linewidth=0.8, linetype="dashed") +
  annotate("text", x=4.3,  y=45, label="4x", colour="#8892A4", size=3.5) +
  annotate("text", x=10.3, y=45, label="10x", colour="#02C39A", size=3.5) +
  annotate("text", x=35,   y=11.5, label="10% threshold", colour="#F59E0B", size=3.5) +
  labs(title="Minimum detectable methylation difference by coverage (80% power)",
       subtitle="Score test, alpha=0.01, baseline methylation=64%",
       x="Reads per cytosine", y="Minimum detectable difference (%)") +
  theme_bw(base_size=12)
ggsave(file.path(OUT_DIR, "coverage_power_curve.pdf"), p2, width=9, height=5)

message("\nDone. Results in: ", OUT_DIR)
