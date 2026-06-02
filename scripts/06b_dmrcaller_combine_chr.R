#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

# Combine per-chromosome DMR results into genome-wide GRanges.
# Run this after all 06_dmrcaller_by_chr.R jobs finish.
# Can also run on partial results to get a preview — missing chromosomes are skipped.

CHROMS <- paste0("chr", c(1:22, "X", "Y"))

CONTRASTS <- c(
  "ASO_VPA_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_CTRL_vs_Scramble_CTRL",
  "ASO_VPA_vs_ASO_CTRL",
  "ASO_VPA_vs_Scramble_VPA"
)

OUT_DIR <- "results/dmr"

for (contrast in CONTRASTS) {
  message("\nCombining: ", contrast)
  chr_dir <- file.path(OUT_DIR, "by_chr", contrast)

  dmr_list <- lapply(CHROMS, function(chr) {
    f <- file.path(chr_dir, paste0("dmr_", contrast, "_", chr, ".rds"))
    if (!file.exists(f)) {
      message("  missing: ", chr)
      return(NULL)
    }
    readRDS(f)
  })

  dmr_list <- Filter(Negate(is.null), dmr_list)
  dmr_list <- Filter(function(x) length(x) > 0, dmr_list)

  if (length(dmr_list) == 0) {
    message("  no results found, skipping")
    next
  }

  all_dmrs <- do.call(c, dmr_list)
  n_total  <- length(all_dmrs)

  # DMRcaller direction note:
  # "gain" = treatment is hypomethylated (proportion1 < proportion2)
  # "loss" = treatment is hypermethylated
  # confirmed by checking proportion1 vs proportion2 in chr13 DMRs
  n_hypo  <- sum(all_dmrs$regionType == "gain",  na.rm=TRUE)
  n_hyper <- sum(all_dmrs$regionType == "loss", na.rm=TRUE)
  message("  total DMRs: ", n_total,
          " (", n_hyper, " hyper / ", n_hypo, " hypo)")
  message("  chr13 DMRs: ",
          sum(as.character(seqnames(all_dmrs)) == "chr13"))
  message("  chromosomes with DMRs: ", length(dmr_list), "/", length(CHROMS))

  saveRDS(all_dmrs, file.path(OUT_DIR, paste0("dmr_", contrast, ".rds")))

  # BED output for IGV — 0-based coordinates, score = -log10(pValue)
  bed_df <- data.frame(
    chr      = as.character(seqnames(all_dmrs)),
    start    = start(all_dmrs) - 1,
    end      = end(all_dmrs),
    name     = contrast,
    score    = round(-log10(all_dmrs$pValue + 1e-300), 2),
    strand   = ".",
    methDiff = round(as.integer(all_dmrs$regionType == "loss") -
                     as.integer(all_dmrs$regionType == "gain"), 4),
    type     = all_dmrs$regionType,
    nCpG     = all_dmrs$cytosinesCount,
    pValue   = all_dmrs$pValue
  )
  write.table(bed_df,
              file.path(OUT_DIR, paste0("dmr_", contrast, ".bed")),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  message("  saved to: ", OUT_DIR)
}
message("\ndone.")
