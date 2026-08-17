#!/usr/bin/env Rscript
# 24_mingap_sensitivity.R
# minGap sensitivity sweep across all chromosomes, four pairwise contrasts.
# Tests minGap 100, 200, 300 (locked), 500, 1000.
# Vectorized: loads and pools each chromosome ONCE per contrast, then
# reuses the pooled data across all 5 minGap values (was reloading per
# minGap value before, causing 5x redundant file I/O).

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/mingap_sensitivity"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

BY_CHR <- "results/alignments/bs/by_chr"
BIN_SIZE <- 300
MIN_DIFF <- 0.20
P_VAL    <- 0.01
MIN_CYT  <- 3
MIN_READS <- 3
MINGAP_VALUES <- c(100, 200, 300, 500, 1000)
KEEP_CHR <- paste0("chr", c(1:22, "X"))

CONDITIONS <- list(
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = paste0("ASO_VPA_",       1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)

CONTRASTS <- list(
  list(name="ASO alone",  a="ASO_CTRL",    b="Scramble_CTRL"),
  list(name="VPA alone",  a="Scramble_VPA",b="Scramble_CTRL"),
  list(name="ASO in VPA", a="ASO_VPA",     b="Scramble_VPA"),
  list(name="VPA in ASO", a="ASO_VPA",     b="ASO_CTRL")
)

load_chr <- function(cond, chr) {
  grs <- lapply(CONDITIONS[[cond]], function(s) {
    f <- file.path(BY_CHR, paste0(s, "_", chr, ".CpG_report.txt.gz"))
    readBismark(f)
  })
  poolMethylationDatasets(GRangesList(grs))
}

results <- list()

for (ct in CONTRASTS) {
  message("\n=== Contrast: ", ct$name, " ===")
  # accumulate counts per minGap across chromosomes
  counts_by_mingap <- setNames(rep(0, length(MINGAP_VALUES)), as.character(MINGAP_VALUES))

  for (chr in KEEP_CHR) {
    message("  ", chr, " - loading once...")
    ma <- tryCatch(load_chr(ct$a, chr), error=function(e) NULL)
    mb <- tryCatch(load_chr(ct$b, chr), error=function(e) NULL)
    if (is.null(ma) || is.null(mb)) {
      message("    skip (load failed)")
      next
    }

    for (mg in MINGAP_VALUES) {
      n <- tryCatch({
        dmrs <- computeDMRs(ma, mb,
                            context="CG", method="bins",
                            binSize=BIN_SIZE,
                            minReadsPerCytosine=MIN_READS,
                            minCytosinesCount=MIN_CYT,
                            minProportionDifference=MIN_DIFF,
                            pValueThreshold=P_VAL,
                            minGap=mg)
        length(dmrs)
      }, error=function(e) { message("    minGap=", mg, " error: ", e$message); 0 })

      counts_by_mingap[as.character(mg)] <- counts_by_mingap[as.character(mg)] + n
      message("    minGap=", mg, ": +", n, " (running total ", counts_by_mingap[as.character(mg)], ")")
    }
    rm(ma, mb); gc(verbose=FALSE)
  }

  for (mg in MINGAP_VALUES) {
    results[[paste(ct$name, mg)]] <- data.frame(
      contrast=ct$name, minGap=mg, n_dmrs=counts_by_mingap[as.character(mg)]
    )
  }

  # checkpoint save after each contrast in case of timeout
  df_partial <- do.call(rbind, results)
  write.csv(df_partial, file.path(OUT, "mingap_sensitivity_partial.csv"), row.names=FALSE)
  message("  Checkpoint saved after contrast: ", ct$name)
}

df <- do.call(rbind, results)
write.csv(df, file.path(OUT, "mingap_sensitivity_full.csv"), row.names=FALSE)
cat("\nFinal results:\n"); print(df)

df$contrast <- factor(df$contrast,
  levels=c("ASO alone","VPA alone","ASO in VPA","VPA in ASO"))

COND_COLS <- c("ASO alone"="#1F3A5F","VPA alone"="#F0A500",
               "ASO in VPA"="#C0392B","VPA in ASO"="#6B4E9E")

p <- ggplot(df, aes(x=minGap, y=n_dmrs, colour=contrast, group=contrast)) +
  geom_line(linewidth=1) +
  geom_point(size=3) +
  geom_vline(xintercept=300, linetype="dashed", colour="grey40") +
  annotate("text", x=310, y=max(df$n_dmrs)*0.95,
           label="Locked\nminGap=300", hjust=0, size=3, colour="grey30") +
  scale_colour_manual(values=COND_COLS, name=NULL) +
  scale_x_continuous(breaks=MINGAP_VALUES) +
  labs(x="minGap (bp)", y="DMR count (genome-wide)") +
  theme_classic(base_size=13) +
  theme(legend.position="top")

ggsave(file.path(OUT, "FigC.3_mingap_sensitivity.pdf"),
       p, width=8, height=5, device=cairo_pdf)
message("Done.")
