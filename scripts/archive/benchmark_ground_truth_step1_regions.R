.libPaths("~/R/library")
library(DMRcaller)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/benchmark_ground_truth"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Loading chr1 ASO_VPA data...")
dat <- readBismarkPool(c(
  file.path(COV_DIR, "ASO_VPA_1_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_2_chr1.CpG_report.txt.gz"),
  file.path(COV_DIR, "ASO_VPA_3_chr1.CpG_report.txt.gz")
))
dat <- dat[mcols(dat)$readsN >= 4]
message("Covered CpGs: ", length(dat))

positions  <- start(dat)
sizes      <- c(100, 500, 1000, 2000)
n_per_size <- 50
spike_regions <- list()

set.seed(42)
for (sz in sizes) {
  count <- 0; attempts <- 0
  while (count < n_per_size && attempts < 100000) {
    attempts  <- attempts + 1
    idx       <- sample(length(positions) - 10, 1)
    start_pos <- positions[idx]
    end_pos   <- start_pos + sz
    cpgs_in   <- sum(positions >= start_pos & positions <= end_pos)
    if (cpgs_in >= 4) {
      count <- count + 1
      spike_regions[[length(spike_regions) + 1]] <- data.frame(
        chr="chr1", start=start_pos, end=end_pos, size=sz, n_cpgs=cpgs_in)
    }
  }
  message("Size ", sz, "bp: ", count, " regions in ", attempts, " attempts")
}

spike_df <- do.call(rbind, spike_regions)
write.csv(spike_df, file.path(OUT_DIR, "spike_in_regions.csv"), row.names=FALSE)
message("Saved ", nrow(spike_df), " regions")
print(table(spike_df$size))
