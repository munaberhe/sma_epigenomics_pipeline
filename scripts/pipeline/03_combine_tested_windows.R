#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({ library(GenomicRanges) })

# Combines per-chromosome tested_windows_<contrast>_<chr>.rds files into
# one genome-wide tested_windows_<contrast>.rds per contrast, matching the
# path that 07_grant_fig1cd_annotation_obsexp.R expects.

TW_DIR  <- "results/dmr/tested_windows"
OUT_DIR <- "results/dmr"

CONTRASTS <- c(
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_CTRL_vs_Scramble_CTRL",
  "ASO_VPA_vs_ASO_CTRL",
  "ASO_VPA_vs_Scramble_VPA"
)
CHROMS <- c(paste0("chr", 1:22), "chrX", "chrY")

for (contrast in CONTRASTS) {
  message("\nCombining: ", contrast)
  files <- file.path(TW_DIR, paste0("tested_windows_", contrast, "_", CHROMS, ".rds"))
  present <- file.exists(files)
  message("  ", sum(present), " / ", length(files), " chromosome files found")
  if (sum(present) == 0) {
    message("  SKIP -- no files for this contrast")
    next
  }
  if (any(!present)) {
    message("  Missing: ", paste(CHROMS[!present], collapse=", "))
  }
  grs <- lapply(files[present], readRDS)
  combined <- do.call(c, grs)
  message("  Total tested windows: ", length(combined))

  out_path <- file.path(OUT_DIR, paste0("tested_windows_", contrast, ".rds"))
  saveRDS(combined, out_path)
  message("  Saved: ", out_path)
}

message("\nDone.")
