#!/usr/bin/env Rscript
# 23_aso_vpa2_sensitivity.R
# Sensitivity analysis: re-run DMR calling excluding ASO_VPA_2.
# ASO_VPA pooled from replicates 1 and 3 only (n=2).
# Compares DMR counts against the original n=3 results.

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(ggrepel)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/dmr_sensitivity_n2"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

BY_CHR_DIR <- "results/alignments/bs/by_chr"
BIN_SIZE <- 300
MIN_DIFF <- 0.20
P_VAL    <- 0.01
MIN_CYT  <- 3
MIN_READS <- 3

KEEP_CHR <- paste0("chr", c(1:22, "X"))

# ASO_VPA_2 excluded
CONDITIONS_N2 <- list(
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = c("ASO_VPA_1", "ASO_VPA_3"),  # ASO_VPA_2 excluded
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)

CONTRASTS <- list(
  ASO_alone  = c("ASO_CTRL",   "Scramble_CTRL"),
  VPA_alone  = c("Scramble_VPA","Scramble_CTRL"),
  ASO_in_VPA = c("ASO_VPA",    "Scramble_VPA"),
  VPA_in_ASO = c("ASO_VPA",    "ASO_CTRL")
)

results_n2 <- list()

for (ct_name in names(CONTRASTS)) {
  out_rds <- file.path(OUT, paste0("dmr_n2_", ct_name, ".rds"))
  if (file.exists(out_rds)) {
    message("Skipping ", ct_name, " - already done")
    results_n2[[ct_name]] <- length(readRDS(out_rds))
    next
  }
  message("\nContrast: ", ct_name)
  cond_a <- CONTRASTS[[ct_name]][1]
  cond_b <- CONTRASTS[[ct_name]][2]
  dmr_list <- list()

  for (chr in KEEP_CHR) {
    message("  ", chr)
    tryCatch({
      load_chr <- function(cond) {
        samples <- CONDITIONS_N2[[cond]]
        grs <- lapply(samples, function(s) {
          f <- file.path(BY_CHR_DIR, paste0(s, "_", chr, ".CpG_report.txt.gz"))
          if (!file.exists(f)) stop("Missing: ", f)
          readBismark(f)
        })
        poolMethylationDatasets(GRangesList(grs))
      }
      ma <- load_chr(cond_a)
      mb <- load_chr(cond_b)
      dmrs <- computeDMRs(ma, mb,
                          context="CG",
                          method="bins",
                          binSize=BIN_SIZE,
                          minReadsPerCytosine=MIN_READS,
                          minCytosinesCount=MIN_CYT,
                          minProportionDifference=MIN_DIFF,
                          pValueThreshold=P_VAL)
      if (length(dmrs) > 0) dmr_list[[chr]] <- dmrs
    }, error=function(e) message("  Error: ", e$message))
  }

  all_dmrs <- if (length(dmr_list) > 0) do.call(c, unname(dmr_list)) else GRanges()
  n <- length(all_dmrs)
  message("  Total DMRs: ", n)
  saveRDS(all_dmrs, file.path(OUT, paste0("dmr_n2_", ct_name, ".rds")))
  results_n2[[ct_name]] <- n
}

# comparison
orig <- c(ASO_alone=3423, VPA_alone=598485, ASO_in_VPA=23669, VPA_in_ASO=664202)
labels <- c(ASO_alone="ASO alone", VPA_alone="VPA alone",
            ASO_in_VPA="ASO in VPA", VPA_in_ASO="VPA in ASO")

df <- data.frame(
  contrast = names(orig),
  label    = labels[names(orig)],
  original = as.numeric(orig),
  n2       = as.numeric(unlist(results_n2[names(orig)]))
)
df$pct_change <- round(100*(df$n2 - df$original)/df$original, 1)
write.csv(df, file.path(OUT, "sensitivity_n2_counts.csv"), row.names=FALSE)
cat("\nResults:\n"); print(df)

p <- ggplot(df, aes(x=original, y=n2, label=label)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") +
  geom_point(colour="#1F3A5F", size=4) +
  geom_text_repel(size=3.5, colour="grey20") +
  scale_x_log10(labels=scales::comma) +
  scale_y_log10(labels=scales::comma) +
  labs(x="DMR count (n=3, all replicates)",
       y="DMR count (ASO_VPA_2 excluded, n=2)",
       caption="Dashed line: identity. Points near line indicate robustness to ASO_VPA_2 exclusion.") +
  theme_classic(base_size=13)

ggsave(file.path(OUT, "FigF.6_sensitivity_n2.pdf"),
       p, width=7, height=6, device=cairo_pdf)
message("Done. Output: ", OUT)
