# ── DMRcaller Differential Methylation Analysis ───────────────────────────
# ASO1+VPA vs ASO1 — SMA Epigenomics Project

library(DMRcaller)
library(BSseq)
library(ggplot2)
library(RColorBrewer)

# ── 1. Setup ──────────────────────────────────────────────────────────────
setwd("~/sma_epigenomics_pipeline")

# Update these with your actual file paths when data arrives
coverage_files <- list.files("results/alignments/bs",
                             pattern = ".bismark.cov.gz",
                             full.names = TRUE)
context  <- "CpG"
min_cov  <- 10

# ── 2. Load bismark coverage files ────────────────────────────────────────
n            <- length(coverage_files)
sample_names <- gsub("_bismark.bismark.cov.gz", "", basename(coverage_files))
conditions   <- c(rep("ASO1", n/2), rep("ASO1_VPA", n/2))

bs <- read.bismark(
  files       = coverage_files,
  sampleNames = sample_names,
  verbose     = TRUE
)
bs

# ── 3. Filter by coverage ─────────────────────────────────────────────────
bs.filtered <- bs[which(rowSums(getCoverage(bs) >= min_cov) == n), ]
cat("CpG sites after filtering:", nrow(bs.filtered), "\n")

# ── 4. Pool samples by condition ──────────────────────────────────────────
pool.ASO1 <- computeMethylationDataPool(bs.filtered[, conditions == "ASO1"])
pool.VPA  <- computeMethylationDataPool(bs.filtered[, conditions == "ASO1_VPA"])

# ── 5. Call DMRs ──────────────────────────────────────────────────────────
# Update method parameter after supervisor meeting
dmrs <- computeDMRs(
  methylationData1        = pool.ASO1,
  methylationData2        = pool.VPA,
  regions                 = NULL,
  context                 = context,
  method                  = "neighbourhood",
  kernelFunction          = "triangular",
  w                       = 3,
  n                       = 3,
  test                    = "score",
  pValueThreshold         = 0.05,
  minCytosines            = 5,
  minProportionDifference = 0.1,
  verbose                 = TRUE
)

dmrs_df <- as.data.frame(dmrs)
cat("DMRs identified:", nrow(dmrs_df), "\n")
head(dmrs_df)

write.csv(dmrs_df, "results/differential/dmrs.csv", row.names = FALSE)

# ── 6. Methylation difference histogram ───────────────────────────────────
ggplot(dmrs_df, aes(x = methylationDifference)) +
  geom_histogram(bins = 50, fill = "#00897B", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  theme_bw() +
  labs(title    = "Distribution of Methylation Differences",
       subtitle = "ASO1+VPA vs ASO1",
       x        = "Methylation Difference",
       y        = "Count")

ggsave("results/figures/dmr_difference_histogram.png", width = 10, height = 7, dpi = 150)

# ── 7. DMR size distribution ──────────────────────────────────────────────
dmrs_df$width <- dmrs_df$end - dmrs_df$start

ggplot(dmrs_df, aes(x = width)) +
  geom_histogram(bins = 50, fill = "#0D1B2A", color = "white") +
  theme_bw() +
  labs(title = "DMR Size Distribution",
       x     = "DMR Width (bp)",
       y     = "Count")

ggsave("results/figures/dmr_size_distribution.png", width = 10, height = 7, dpi = 150)

# ── 8. DMRs per chromosome ────────────────────────────────────────────────
chr_counts <- as.data.frame(table(dmrs_df$seqnames))
colnames(chr_counts) <- c("chromosome", "count")

ggplot(chr_counts, aes(x = chromosome, y = count)) +
  geom_bar(stat = "identity", fill = "#00897B") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of DMRs per Chromosome",
       x     = "Chromosome",
       y     = "Number of DMRs")

ggsave("results/figures/dmr_per_chromosome.png", width = 12, height = 7, dpi = 150)

# ── 9. Genome-wide methylation per sample ─────────────────────────────────
meth_levels <- data.frame(
  sample    = sample_names,
  condition = conditions,
  mean_meth = colMeans(getMeth(bs.filtered, type = "raw"), na.rm = TRUE)
)

ggplot(meth_levels, aes(x = sample, y = mean_meth, fill = condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("ASO1" = "#0D1B2A", "ASO1_VPA" = "#00897B")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Genome-Wide Mean Methylation per Sample",
       x     = "Sample",
       y     = "Mean Methylation Level")

ggsave("results/figures/genome_wide_methylation.png", width = 10, height = 7, dpi = 150)

cat("DMRcaller analysis complete. All plots saved to results/figures/\n")