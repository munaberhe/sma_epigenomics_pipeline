#!/usr/bin/env Rscript
# 21_enhancer_enrichment_pairwise.R
# Genome-wide DMR enrichment at H9 predicted enhancers, four pairwise contrasts.
# Method: regioneR permutation test (1000 permutations, chr-matched random regions).
# Caveat: H9 is an ESC line; HEK293T is not ESC-derived. Enrichment reflects
# overlap with regions that were enhancer-active in stem cells.
# Input: results/dmr/dmr_<contrast>.rds (locked params: binSize=300, minDiff=0.20, p<0.01)
# Output: results/enhancer_pairwise/

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(regioneR)
  library(ggplot2)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/enhancer_pairwise"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

N_PERM   <- 1000
KEEP_CHR <- paste0("chr", c(1:22, "X"))
ENH_BED  <- "data/reference/H9_predicted_non_promoter_non_fragments.bed.gz"

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",   label="ASO alone"),
  list(name="Scramble_VPA_vs_Scramble_CTRL", label="VPA alone"),
  list(name="ASO_VPA_vs_Scramble_VPA",     label="ASO in VPA"),
  list(name="ASO_VPA_vs_ASO_CTRL",         label="VPA in ASO")
)

message("Loading H9 enhancers...")
bed <- read.table(pipe(paste0("zcat ", ENH_BED)),
                  header=TRUE, sep="\t", stringsAsFactors=FALSE)
bed <- bed[bed$seqnames %in% KEEP_CHR, ]
enh_gr <- GRanges(bed$seqnames, IRanges(bed$start, bed$end))
message("  Enhancers: ", length(enh_gr))

results <- list()

for (ct in CONTRASTS) {
  message("Running: ", ct$label)
  dmrs <- readRDS(paste0("results/dmr/dmr_", ct$name, ".rds"))
  dmrs <- dmrs[seqnames(dmrs) %in% KEEP_CHR]
  # subsample all contrasts to 50,000 DMRs for uniform computational treatment
  MAX_DMRS <- 50000
  n_orig <- length(dmrs)
  if (length(dmrs) > MAX_DMRS) {
    set.seed(42)
    dmrs <- dmrs[sample(length(dmrs), MAX_DMRS)]
  }
  message("  DMRs: ", length(dmrs), " (original: ", n_orig, ")")

  pt <- tryCatch(
    permTest(A=dmrs, B=enh_gr,
             randomize.function=randomizeRegions,
             evaluate.function=numOverlaps,
             genome="hg38",
             ntimes=N_PERM,
             verbose=FALSE),
    error=function(e) { message("  Error: ", e$message); NULL }
  )

  if (!is.null(pt)) {
    obs   <- pt$numOverlaps$observed
    nullm <- mean(pt$numOverlaps$permuted)
    pval  <- pt$numOverlaps$pval
    zscore <- pt$numOverlaps$zscore
    fold  <- obs / nullm
    message(sprintf("  Obs=%d Null=%.1f Fold=%.2f p=%.3f z=%.2f",
                    obs, nullm, fold, pval, zscore))
    results[[ct$label]] <- data.frame(
      contrast=ct$label, observed=obs, null_mean=nullm,
      fold_enrichment=fold, pvalue=pval, zscore=zscore
    )
  }
}

df <- do.call(rbind, results)
write.csv(df, file.path(OUT, "enhancer_enrichment_pairwise.csv"), row.names=FALSE)
cat("\nResults:\n")
print(df)

# plot
df$contrast <- factor(df$contrast, levels=sapply(CONTRASTS, function(x) x$label))
df$sig <- ifelse(df$pvalue < 0.05, "*", "")

p <- ggplot(df, aes(x=contrast, y=fold_enrichment, fill=contrast)) +
  geom_col(width=0.6) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey40") +
  geom_text(aes(label=sig), vjust=-0.5, size=6) +
  scale_fill_manual(values=c(
    "ASO alone"    = "#1F3A5F",
    "VPA alone"    = "#F0A500",
    "ASO in VPA"   = "#C0392B",
    "VPA in ASO"   = "#C0392B"
  )) +
  labs(x=NULL, y="Fold enrichment over random",
       caption="H9 predicted non-promoter enhancers. 1000 permutations, chr-matched null.") +
  theme_classic(base_size=13) +
  theme(legend.position="none",
        axis.text.x=element_text(size=11))

ggsave(file.path(OUT, "enhancer_enrichment_pairwise.pdf"),
       p, width=8, height=5, device=cairo_pdf)
message("Done. Output: ", OUT)
