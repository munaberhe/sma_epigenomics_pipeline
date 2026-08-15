.libPaths("~/R/library")
# benchmark_ground_truth_step2_spikein.R
# Step 2: Create spiked datasets for ground truth benchmark
# Copies real data and artificially reduces methylation in treatment
# at spike-in regions. Real data is NEVER modified.
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/benchmark_ground_truth"

# Load spike-in regions
spike_df <- read.csv(file.path(OUT_DIR, "spike_in_regions.csv"))
message("Loaded ", nrow(spike_df), " spike-in regions")

spike_regions <- GRanges(
  seqnames = spike_df$chr,
  ranges   = IRanges(spike_df$start, spike_df$end)
)

# Load real chr1 data — copies only, never modify originals
message("Loading chr1 data...")
aso_vpa <- readBismarkPool(c(
  file.path(COV_DIR, "ASO_VPA_1_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_2_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_3_chr1.CpG_report.txt.gz")
))
aso_ctrl <- readBismarkPool(c(
  file.path(COV_DIR, "ASO_CTRL_1_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_CTRL_2_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_CTRL_3_chr1.CpG_report.txt.gz")
))
message("ASO_VPA CpGs:  ", length(aso_vpa))
message("ASO_CTRL CpGs: ", length(aso_ctrl))

# Function to spike in methylation differences
# Reduces readsM in treatment at spike regions by delta_meth
# Control data is never touched
spike_in <- function(treatment, regions, delta_meth=0.2, seed=42) {
  set.seed(seed)
  treat_spiked <- treatment  # copy — original unchanged
  hits <- findOverlaps(treat_spiked, regions)
  spike_idx <- unique(queryHits(hits))
  message("  CpGs in spike regions: ", length(spike_idx))

  # Reduce methylation proportionally
  original_M <- mcols(treat_spiked)$readsM[spike_idx]
  original_N <- mcols(treat_spiked)$readsN[spike_idx]
  new_M <- round(original_M * (1 - delta_meth))
  # ensure readsM never exceeds readsN and never goes below 0
  new_M <- pmax(0, pmin(new_M, original_N))
  mcols(treat_spiked)$readsM[spike_idx] <- new_M

  message("  Mean methylation before: ",
          round(sum(original_M)/sum(original_N)*100, 1), "%")
  message("  Mean methylation after:  ",
          round(sum(new_M)/sum(original_N)*100, 1), "%")
  treat_spiked
}

# Create spiked datasets at four delta_meth levels
delta_meths <- c(0.1, 0.2, 0.3, 0.4)

for (dm in delta_meths) {
  message("\n--- Spiking at delta_meth=", dm, " ---")
  aso_vpa_spiked <- spike_in(aso_vpa, spike_regions, delta_meth=dm, seed=42)

  # Save spiked treatment and original control as RDS
  out_treat <- file.path(OUT_DIR,
    paste0("aso_vpa_spiked_dm", gsub("\\.", "p", dm), ".rds"))
  out_ctrl  <- file.path(OUT_DIR, "aso_ctrl_original.rds")

  saveRDS(aso_vpa_spiked, out_treat)
  message("  Saved spiked treatment: ", out_treat)
}

# Save original control once
saveRDS(aso_ctrl, file.path(OUT_DIR, "aso_ctrl_original.rds"))
message("\nSaved original control: ", OUT_DIR, "/aso_ctrl_original.rds")
message("\nStep 2 complete. Spiked datasets saved to: ", OUT_DIR)
message("Real data untouched — originals in: ", COV_DIR)
