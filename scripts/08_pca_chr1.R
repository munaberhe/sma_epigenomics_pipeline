#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
})

# PCA of per-replicate CpG methylation using chr1 as representative chromosome.
# Each dot = one biological replicate, coloured by condition.
# Also generates per-sample methylation violin plots.

BY_CHR  <- "results/alignments/bs/by_chr"
OUT     <- "results/dmr_qc"
MIN_COV <- 5  # minimum coverage required in all 12 samples
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

SAMPLES <- c(
  "ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3",
  "ASO_VPA_1",  "ASO_VPA_2",  "ASO_VPA_3",
  "Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
  "Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3"
)
CONDITION <- c(rep("ASO_CTRL",3), rep("ASO_VPA",3),
               rep("Scramble_CTRL",3), rep("Scramble_VPA",3))

GROUP_COLS <- c(
  ASO_CTRL      = "#2E9B6F",
  ASO_VPA       = "#D94F3D",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#1D6FA4"
)

# read chr1 CpG reports for all 12 samples
message("reading chr1 for all 12 samples...")
meth_list <- list()
for (s in SAMPLES) {
  f  <- file.path(BY_CHR, paste0(s, "_chr1.CpG_report.txt.gz"))
  gr <- readBismark(f)
  meth_list[[s]] <- gr[gr$context == "CG"]
  message("  ", s)
}

# find CpGs covered at >= MIN_COV in all 12 samples
message("finding common CpGs (>= ", MIN_COV, "x in all samples)...")
common_pos <- Reduce(intersect, lapply(meth_list, function(gr) {
  gr <- gr[gr$readsN >= MIN_COV]
  paste0(as.character(seqnames(gr)), ":", start(gr))
}))
message("  common CpGs: ", length(common_pos))

# build methylation matrix: samples x CpGs
message("building matrix...")
meth_mat <- do.call(cbind, lapply(SAMPLES, function(s) {
  gr  <- meth_list[[s]]
  pos <- paste0(as.character(seqnames(gr)), ":", start(gr))
  idx <- match(common_pos, pos)
  gr$readsM[idx] / gr$readsN[idx]
}))
colnames(meth_mat) <- SAMPLES
rownames(meth_mat) <- common_pos
meth_mat <- meth_mat[complete.cases(meth_mat), ]
message("  matrix: ", nrow(meth_mat), " CpGs x ", ncol(meth_mat), " samples")

# subsample to 100k CpGs for speed
if (nrow(meth_mat) > 100000) {
  set.seed(42)
  meth_mat <- meth_mat[sample(nrow(meth_mat), 100000), ]
  message("  subsampled to 100k CpGs")
}

# PCA — transpose so samples are rows
message("running PCA...")
pca    <- prcomp(t(meth_mat), scale.=FALSE, center=TRUE)
pca_df <- data.frame(
  PC1       = pca$x[,1],
  PC2       = pca$x[,2],
  sample    = SAMPLES,
  condition = CONDITION,
  replicate = rep(c("1","2","3"), 4)
)
var_exp <- round(summary(pca)$importance[2,1:2]*100, 1)
message("PC1: ", var_exp[1], "% | PC2: ", var_exp[2], "%")

p <- ggplot(pca_df, aes(PC1, PC2, colour=condition,
                         shape=replicate, label=sample)) +
  geom_point(size=5, alpha=0.9) +
  geom_text(vjust=-0.9, size=2.8, colour="grey30") +
  stat_ellipse(aes(group=condition), type="norm",
               linetype="dashed", linewidth=0.4, alpha=0.6) +
  scale_colour_manual(values=GROUP_COLS, name="Condition") +
  scale_shape_manual(values=c("1"=16,"2"=17,"3"=15), name="Replicate") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold"), legend.position="right") +
  labs(
    title    = "PCA of per-replicate CpG methylation (chr1)",
    subtitle = paste0("n=12 | ", format(nrow(meth_mat), big.mark=","),
                      " CpGs | min coverage ", MIN_COV, "x"),
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  )
ggsave(file.path(OUT, "sample_PCA_12samples_chr1.pdf"), p, width=9, height=7)
ggsave(file.path(OUT, "sample_PCA_12samples_chr1.png"), p, width=9, height=7, dpi=200)
message("saved PCA plot")

# per-sample methylation violin
message("generating violin plot...")
violin_df <- data.frame(
  methylation = as.vector(meth_mat),
  sample      = rep(colnames(meth_mat), each=nrow(meth_mat))
)
violin_df$condition <- CONDITION[match(violin_df$sample, SAMPLES)]
violin_df$sample    <- factor(violin_df$sample, levels=SAMPLES)
violin_df$condition <- factor(violin_df$condition,
  levels=c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA"))

p_violin <- ggplot(violin_df, aes(x=sample, y=methylation, fill=condition)) +
  geom_violin(trim=FALSE, alpha=0.8, linewidth=0.3) +
  geom_boxplot(width=0.07, fill="white",
               outlier.size=0.2, outlier.alpha=0.2) +
  scale_fill_manual(values=GROUP_COLS) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(title="Per-sample CpG methylation distribution (chr1)",
       subtitle=paste0(format(nrow(meth_mat), big.mark=","),
                       " CpGs at >=5x | VPA samples show global hypomethylation"),
       x=NULL, y="CpG methylation proportion") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        plot.title=element_text(face="bold"), legend.position="top")
ggsave(file.path(OUT, "per_sample_methylation_violin_chr1.pdf"),
       p_violin, width=14, height=6)
ggsave(file.path(OUT, "per_sample_methylation_violin_chr1.png"),
       p_violin, width=14, height=6, dpi=200)
message("saved violin plot")
message("done.")
