# genome-wide DMR vs H9 enhancer enrichment.
# question: do ASO/VPA DMRs land at H9 non-promoter enhancers more (or less)
# often than expected by chance?
# method: count observed DMR-enhancer overlaps per contrast; compare against
# a chr/width-matched random region null (regioneR permutation test).
# caveat: H9 is an ESC line; samples are not ESCs. positive enrichment means
# DMRs hit regions that were enhancer-active in stem cells (often poised or
# silenced in differentiated cells). null does not exclude lineage-specific
# enhancer effects from cell-type-matched maps.

.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(regioneR)
  library(ggplot2)
})
source("scripts/00_sma_palette.R")
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/genomewide_enhancer"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ---- params ----
N_PERM       <- 1000          # permutations per contrast
GENOME       <- "hg38"
ENHANCER_BED <- "data/reference/H9_predicted_non_promoter_non_fragments.bed.gz"
CONTRASTS <- c("ASO_CTRL_vs_Scramble_CTRL","ASO_VPA_vs_Scramble_CTRL",
               "ASO_VPA_vs_ASO_CTRL","ASO_VPA_vs_Scramble_VPA",
               "Scramble_VPA_vs_Scramble_CTRL")
# canonical chromosomes only - permTest works best on standard chroms
KEEP_CHR <- paste0("chr", c(1:22, "X"))

# ---- load enhancers ----
message("loading H9 enhancers...")
bed <- read.table(pipe(paste0("zcat ", ENHANCER_BED)),
                  header=TRUE, sep="\t", stringsAsFactors=FALSE)
bed <- bed[bed$seqnames %in% KEEP_CHR, ]
enh_gr <- GRanges(bed$seqnames, IRanges(bed$start, bed$end), type=bed$type)
message("  enhancers loaded: ", length(enh_gr),
        " (on ", length(unique(seqnames(enh_gr))), " chromosomes)")

# pre-compute total enhancer span and genome span for the expected baseline.
# enhancer span / genome span = expected fraction of DMRs overlapping by chance
# (if DMRs were placed uniformly at random with no chr bias).
total_enh_span <- sum(as.numeric(width(enh_gr)))
message("  total enhancer span (bp): ", format(total_enh_span, big.mark=","))

# ---- per-contrast: observed overlaps + permutation null ----
results <- list()

for (ct in CONTRASTS) {
  message("\n=== ", ct, " ===")
  dmr <- readRDS(paste0("results/dmr/dmr_", ct, ".rds"))

  # keep canonical chr only, dedupe by coordinate key
  dmr <- dmr[as.character(seqnames(dmr)) %in% KEEP_CHR]
  key <- paste(seqnames(dmr), start(dmr), end(dmr), sep="_")
  dmr <- dmr[!duplicated(key)]
  message("  DMRs (canonical chr, dedup): ", length(dmr))

  if (length(dmr) == 0) {
    message("  no DMRs, skipping")
    next
  }

  # observed overlap count
  obs_hits <- length(subsetByOverlaps(dmr, enh_gr))
  obs_frac <- obs_hits / length(dmr)
  message("  observed: ", obs_hits, " / ", length(dmr),
          " DMRs overlap an enhancer (", round(100*obs_frac, 2), "%)")

  # permutation test - shuffle DMRs across the genome preserving chr and width
  message("  running ", N_PERM, " permutations...")
  set.seed(42)
  pt <- permTest(
    A                   = dmr,
    B                   = enh_gr,
    ntimes              = N_PERM,
    randomize.function  = randomizeRegions,
    evaluate.function   = numOverlaps,
    genome              = GENOME,
    allow.overlaps      = TRUE,
    per.chromosome      = TRUE,    # preserve chr distribution of DMRs
    mask                = NULL,
    count.once          = TRUE,    # count each DMR once even if it hits multiple enhancers
    verbose             = FALSE
  )

  pv  <- pt$numOverlaps$pval
  z   <- pt$numOverlaps$zscore
  exp_mean <- mean(pt$numOverlaps$permuted)
  fold     <- if (exp_mean > 0) obs_hits / exp_mean else NA_real_

  message(sprintf("  permuted mean: %.1f | observed: %d | fold: %.2fx | z=%.2f | p=%.3g",
                  exp_mean, obs_hits, fold, z, pv))

  # save the permutation distribution for plotting
  results[[ct]] <- list(
    contrast    = ct,
    n_dmr       = length(dmr),
    n_observed  = obs_hits,
    n_expected  = exp_mean,
    fold        = fold,
    zscore      = z,
    pvalue      = pv,
    permuted    = pt$numOverlaps$permuted
  )

  # per-contrast plot of the null distribution + observed line
  perm_df <- data.frame(n = pt$numOverlaps$permuted)
  p_perm <- ggplot(perm_df, aes(x=n)) +
    geom_histogram(bins=40, fill=SMA_PALETTE$null_fill, colour=SMA_PALETTE$null_outline) +
    geom_vline(xintercept=obs_hits, colour=SMA_PALETTE$observed, linewidth=1.2) +
    annotate("text", x=obs_hits, y=Inf,
             label=sprintf("observed = %d\nfold = %.2fx\np = %.3g",
                           obs_hits, fold, pv),
             hjust=ifelse(obs_hits > median(perm_df$n), 1.05, -0.05),
             vjust=1.3, size=3.2, colour="#B2182B") +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold")) +
    labs(title=paste0("DMR-enhancer overlap permutation null: ", ct),
         subtitle=sprintf(
           "%d DMRs vs %d H9 enhancers | %d permutations | per-chr randomisation",
           length(dmr), length(enh_gr), N_PERM),
         x="number of DMRs overlapping an enhancer (per permutation)",
         y="permutation count")
  ggsave(file.path(OUT_DIR, paste0("perm_null_", ct, ".pdf")),
         p_perm, width=8, height=4.5, device=cairo_pdf)
  saveRDS(results[[ct]], file.path(OUT_DIR, paste0("checkpoint_", ct, ".rds")))
  message("  checkpoint saved: ", ct)
}

# ---- summary table ----
if (length(results) == 0) {
  stop("no contrasts produced results")
}

sum_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    contrast    = r$contrast,
    n_dmr       = r$n_dmr,
    n_observed  = r$n_observed,
    n_expected  = round(r$n_expected, 1),
    fold        = round(r$fold, 3),
    zscore      = round(r$zscore, 2),
    pvalue      = r$pvalue,
    direction   = ifelse(r$fold > 1, "enriched", "depleted"),
    stringsAsFactors = FALSE
  )
}))
sum_df$padj <- p.adjust(sum_df$pvalue, method="BH")

print(sum_df, row.names=FALSE)
write.csv(sum_df,
          file.path(OUT_DIR, "genomewide_enhancer_enrichment_summary.csv"),
          row.names=FALSE)
message("\nsaved: genomewide_enhancer_enrichment_summary.csv")

# ---- combined fold-enrichment plot ----
sum_df$contrast <- factor(sum_df$contrast, levels=CONTRASTS)
sum_df$sig <- ifelse(sum_df$padj < 0.05, "p.adj < 0.05", "ns")

p_fold <- ggplot(sum_df, aes(x=contrast, y=fold, fill=sig)) +
  geom_col(width=0.7, colour="grey20") +
  geom_hline(yintercept=1, linetype="dashed", colour="grey40") +
  geom_text(aes(label=sprintf("%.2fx\np=%.2g", fold, padj)),
            vjust=ifelse(sum_df$fold > 1, -0.3, 1.3),
            size=3.0) +
  scale_fill_manual(values=c("ns"="grey80", "p.adj < 0.05"="#B2182B")) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        plot.title=element_text(face="bold"),
        legend.position="top") +
  labs(title="DMR-enhancer overlap, fold vs chr-matched random",
       subtitle=sprintf(
         "H9 ESC predicted non-promoter enhancers (n=%d) | %d permutations | BH FDR",
         length(enh_gr), N_PERM),
       x=NULL, y="fold enrichment (observed / expected)")

ggsave(file.path(OUT_DIR, "genomewide_enhancer_fold_enrichment.pdf"),
       p_fold, width=9, height=5, device=cairo_pdf)
message("saved: genomewide_enhancer_fold_enrichment.pdf")
message("\ndone. outputs in: ", OUT_DIR)