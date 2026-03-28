# DMRcaller differential methylation analysis
# SMA epigenomics project - Muna Berhe
# Comparing ASO1+VPA vs ASO1 treatment conditions

library(DMRcaller)
library(BSseq)
library(ggplot2)
library(RColorBrewer)

setwd("~/sma_epigenomics_pipeline")

# parameters - method and thresholds to be confirmed with supervisor
context <- "CpG"
min_cov <- 10

# find bismark coverage files
coverage_files <- list.files("results/alignments/bs",
                             pattern = ".bismark.cov.gz",
                             full.names = TRUE)

n            <- length(coverage_files)
sample_names <- gsub("_bismark.bismark.cov.gz", "", basename(coverage_files))
conditions   <- c(rep("ASO1", n/2), rep("ASO1_VPA", n/2))

# load coverage files into BSseq object
bs <- read.bismark(
  files       = coverage_files,
  sampleNames = sample_names,
  verbose     = TRUE
)
bs

# filter sites with insufficient coverage across all samples
bs.filtered <- bs[which(rowSums(getCoverage(bs) >= min_cov) == n), ]
cat("CpG sites after coverage filtering:", nrow(bs.filtered), "\n")

# pool replicates by condition for DMR calling
pool.ASO1 <- computeMethylationDataPool(bs.filtered[, conditions == "ASO1"])
pool.VPA  <- computeMethylationDataPool(bs.filtered[, conditions == "ASO1_VPA"])

# call DMRs using neighbourhood method
# w and n control the kernel width - see DMRcaller vignette for details
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

# distribution of methylation differences
ggplot(dmrs_df, aes(x = methylationDifference)) +
  geom_histogram(bins = 50, fill = "#00897B", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  theme_bw() +
  labs(title = "Distribution of methylation differences",
       subtitle = "ASO1+VPA vs ASO1",
       x = "Methylation difference",
       y = "Count")

ggsave("results/figures/dmr_difference_histogram.png", width = 10, height = 7, dpi = 150)

# DMR size distribution
dmrs_df$width <- dmrs_df$end - dmrs_df$start

ggplot(dmrs_df, aes(x = width)) +
  geom_histogram(bins = 50, fill = "#0D1B2A", color = "white") +
  theme_bw() +
  labs(title = "DMR size distribution",
       x = "DMR width (bp)",
       y = "Count")

ggsave("results/figures/dmr_size_distribution.png", width = 10, height = 7, dpi = 150)

# DMRs per chromosome - check for any unexpected chromosomal enrichment
chr_counts <- as.data.frame(table(dmrs_df$seqnames))
colnames(chr_counts) <- c("chromosome", "count")

ggplot(chr_counts, aes(x = chromosome, y = count)) +
  geom_bar(stat = "identity", fill = "#00897B") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "DMRs per chromosome",
       x = "Chromosome",
       y = "Number of DMRs")

ggsave("results/figures/dmr_per_chromosome.png", width = 12, height = 7, dpi = 150)

# genome-wide mean methylation per sample - check consistency within groups
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
  labs(title = "Genome-wide mean methylation per sample",
       x = "Sample",
       y = "Mean methylation level")

ggsave("results/figures/genome_wide_methylation.png", width = 10, height = 7, dpi = 150)
