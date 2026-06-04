#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(SummarizedExperiment)
})

setwd("/data/home/bt25018/sma_epigenomics_pipeline")

# Combined multi-panel thesis figures from saved RDS outputs.
# All figures sized for thesis submission (A4/letter).
# Loads pre-computed data — does not re-run any analysis.

OUT <- "results/thesis_figures"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

COLS <- c(
  ASO_CTRL      = "#1B4F8A",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#6B7280",
  Hypo          = "#4393C3",
  Hyper         = "#D6604D"
)

CHROMS <- paste0("chr", c(1:22, "X", "Y"))

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",      label="ASO alone vs Scramble CTRL"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",  label="VPA alone vs Scramble CTRL"),
  list(name="ASO_VPA_vs_ASO_CTRL",            label="ASO+VPA vs ASO alone")
)

message("loading DMR objects...")
dmrs <- lapply(CONTRASTS, function(ct) {
  d <- readRDS(paste0("results/dmr/dmr_", ct$name, ".rds"))
  d[mcols(d)$cytosinesCount >= 6]
})
names(dmrs) <- sapply(CONTRASTS, `[[`, "name")

# Figure 1: methylation difference histograms
message("figure 1: methylation diff histograms...")
make_hist <- function(dmr_obj, label, tag) {
  df <- data.frame(
    methDiff  = mcols(dmr_obj)$proportion1 - mcols(dmr_obj)$proportion2,
    direction = ifelse(mcols(dmr_obj)$regionType=="gain", "Hypo", "Hyper")
  )
  ggplot(df, aes(x=methDiff, fill=direction)) +
    geom_histogram(bins=60, colour="white", linewidth=0.1) +
    geom_vline(xintercept=0, linewidth=0.5, linetype="dashed", colour="grey20") +
    geom_vline(xintercept=c(-0.2,0.2), linewidth=0.3, linetype="dotted", colour="grey50") +
    scale_fill_manual(values=c(Hypo=COLS["Hypo"], Hyper=COLS["Hyper"]), name=NULL) +
    scale_x_continuous(limits=c(-1,1), breaks=seq(-1,1,0.5),
                       labels=scales::percent_format(accuracy=1)) +
    scale_y_continuous(labels=scales::comma) +
    labs(title=paste0("(", tag, ") ", label),
         x="methylation difference (treatment minus reference)", y="DMRs") +
    theme_bw(base_size=10) +
    theme(plot.title=element_text(face="bold", size=9),
          legend.position="top", panel.grid.minor=element_blank())
}
fig1 <- (make_hist(dmrs[[1]], CONTRASTS[[1]]$label, "A") |
          make_hist(dmrs[[2]], CONTRASTS[[2]]$label, "B") |
          make_hist(dmrs[[3]], CONTRASTS[[3]]$label, "C")) +
  plot_annotation(title="Figure 1 — Methylation difference distributions",
                  theme=theme(plot.title=element_text(face="bold", size=11)))
ggsave(file.path(OUT, "Fig1_methylation_diff_histograms.pdf"), fig1, width=15, height=5)
message("saved Fig1")

# Figure 2: per-chromosome bar charts
message("figure 2: per-chromosome bars...")
make_chr_bar <- function(dmr_obj, label, tag) {
  chr_counts <- as.data.frame(table(
    chr       = as.character(seqnames(dmr_obj)),
    direction = ifelse(mcols(dmr_obj)$regionType=="gain", "Hypo", "Hyper")
  ))
  chr_counts$chr <- factor(chr_counts$chr, levels=CHROMS)
  chr_counts <- chr_counts[!is.na(chr_counts$chr), ]
  chr_counts$Freq_signed <- ifelse(chr_counts$direction=="Hypo",
                                   -chr_counts$Freq, chr_counts$Freq)
  ggplot(chr_counts, aes(x=chr, y=Freq_signed, fill=direction)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept=0, linewidth=0.4, colour="grey30") +
    scale_fill_manual(values=c(Hypo=COLS["Hypo"], Hyper=COLS["Hyper"]), name=NULL) +
    scale_y_continuous(labels=function(x) format(abs(x), big.mark=",")) +
    labs(title=paste0("(", tag, ") ", label), x="chromosome", y="DMRs") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
          plot.title=element_text(face="bold", size=9),
          legend.position="top", panel.grid.minor=element_blank())
}
fig2 <- (make_chr_bar(dmrs[[1]], CONTRASTS[[1]]$label, "A") /
          make_chr_bar(dmrs[[2]], CONTRASTS[[2]]$label, "B")) +
  plot_annotation(title="Figure 2 — DMR chromosomal distribution",
                  theme=theme(plot.title=element_text(face="bold", size=11)))
ggsave(file.path(OUT, "Fig2_per_chromosome_bars.pdf"), fig2, width=14, height=10)
message("saved Fig2")

# Figure 3: QC triptych
message("figure 3: QC triptych...")
cor_mat  <- read.csv("results/qc/methylation/methylation_correlation.csv", row.names=1)
cor_melt <- reshape2::melt(as.matrix(cor_mat))
p3a <- ggplot(cor_melt, aes(Var1, Var2, fill=value)) +
  geom_tile(colour="white", linewidth=0.3) +
  geom_text(aes(label=round(value,2)), size=2.2, fontface="bold") +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                       midpoint=0.87, limits=c(0.84,1), name="Pearson r") +
  labs(title="(A) Sample correlation heatmap", x=NULL, y=NULL) +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        axis.text.y=element_text(size=7), plot.title=element_text(face="bold", size=9),
        panel.grid=element_blank())

dmr_aso  <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
dmr_aso  <- dmr_aso[mcols(dmr_aso)$cytosinesCount >= 6]
dmr_vpa  <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
dmr_vpa  <- dmr_vpa[mcols(dmr_vpa)$cytosinesCount >= 6]
dmr_comb <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
dmr_comb <- dmr_comb[mcols(dmr_comb)$cytosinesCount >= 6]

violin_df <- rbind(
  data.frame(methylation=mcols(dmr_aso)$proportion1,  condition="ASO_CTRL"),
  data.frame(methylation=mcols(dmr_aso)$proportion2,  condition="Scramble_CTRL"),
  data.frame(methylation=mcols(dmr_vpa)$proportion1,  condition="Scramble_VPA"),
  data.frame(methylation=mcols(dmr_comb)$proportion1, condition="ASO_VPA")
)
violin_df$condition <- factor(violin_df$condition,
  levels=c("ASO_CTRL","Scramble_CTRL","Scramble_VPA","ASO_VPA"))

p3b <- ggplot(violin_df, aes(x=condition, y=methylation, fill=condition)) +
  geom_violin(trim=FALSE, alpha=0.8) +
  geom_boxplot(width=0.07, fill="white", outlier.size=0.2, outlier.alpha=0.2) +
  scale_fill_manual(values=unname(COLS[c("ASO_CTRL","Scramble_CTRL","Scramble_VPA","ASO_VPA")])) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(title="(B) CpG methylation at DMR loci", x=NULL, y="CpG methylation") +
  theme_bw(base_size=9) +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1, size=8),
        plot.title=element_text(face="bold", size=9))

fig3 <- (p3a | p3b) +
  plot_annotation(title="Figure 3 — Quality control",
                  theme=theme(plot.title=element_text(face="bold", size=11)))
ggsave(file.path(OUT, "Fig3_QC.pdf"), fig3, width=14, height=6)
message("saved Fig3")

# Figure 4: H3K9me2
message("figure 4: H3K9me2...")
plot_df <- readRDS("results/h3k9me2_overlap/h3k9me2_plot_df.rds")
plot_df$group <- factor(plot_df$group, levels=c("Background","ASO DMRs"))
p4a <- ggplot(plot_df, aes(x=group, y=signal, fill=group)) +
  geom_boxplot(outlier.size=0.3, outlier.alpha=0.3) +
  facet_wrap(~set, scales="free_y") +
  scale_fill_manual(values=c("Background"="#cccccc","ASO DMRs"="#2E9B6F"), name=NULL) +
  scale_y_continuous(limits=c(0, quantile(plot_df$signal, 0.99, na.rm=TRUE))) +
  labs(title="(A) H3K9me2 at ASO DMRs vs background", x=NULL, y="mean H3K9me2 signal") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=9), legend.position="none",
        strip.text=element_text(face="bold"))

smn2_df <- read.table("results/h3k9me2_overlap/smn2_h3k9me2_signal.tsv",
                      header=TRUE, sep="\t")
smn2_df$treatment <- ifelse(grepl("CTR", smn2_df$condition), "CTR", "ASO")
p4b <- ggplot(smn2_df, aes(x=condition, y=signal, fill=treatment)) +
  geom_col() +
  scale_fill_manual(values=c(CTR="#1D6FA4", ASO="#2E9B6F")) +
  labs(title="(B) H3K9me2 at SMN2 locus", x=NULL, y="mean H3K9me2 signal") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=9))

fig4 <- (p4a | p4b) +
  plot_annotation(title="Figure 4 — H3K9me2 enrichment analysis",
                  theme=theme(plot.title=element_text(face="bold", size=11)))
ggsave(file.path(OUT, "Fig4_H3K9me2.pdf"), fig4, width=12, height=5)
message("saved Fig4")

# Figure 5: negative results
message("figure 5: negative results...")
se <- readRDS("results/tf_motif/motif_enrichment_results.rds")
motif_df <- data.frame(
  motif       = rownames(se),
  log2enr     = as.numeric(assay(se, "log2enr")[,1]),
  neglog10padj = as.numeric(assay(se, "negLog10Padj")[,1])
)
motif_df <- motif_df[complete.cases(motif_df), ]
motif_df$padj <- 10^(-motif_df$neglog10padj)
top5 <- head(motif_df[order(motif_df$padj), ], 5)

p5a <- ggplot(motif_df, aes(x=log2enr, y=neglog10padj)) +
  geom_point(aes(colour=padj<0.05), alpha=0.5, size=1) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="red", linewidth=0.5) +
  scale_colour_manual(values=c("FALSE"="grey70","TRUE"="red"), guide="none") +
  geom_text(data=top5, aes(label=gsub("_.*","",motif)),
            size=2.5, vjust=-0.5, colour="grey30") +
  labs(title="(A) TF motif enrichment — negative result",
       subtitle=NULL,
       x="log2 enrichment", y="-log10(p.adj)") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=9))

sj_df      <- read.csv("results/splice_junction/splice_junction_distances.csv")
sj_summary <- read.csv("results/splice_junction/splice_junction_summary.csv")
p5b <- ggplot(sj_df, aes(x=distance, fill=type, colour=type)) +
  geom_density(alpha=0.4, linewidth=0.7) +
  scale_x_log10(labels=scales::comma) +
  scale_fill_manual(values=c("ASO DMRs"="#2E9B6F","Background"="grey70"), name=NULL) +
  scale_colour_manual(values=c("ASO DMRs"="#2E9B6F","Background"="grey50"), name=NULL) +
  geom_vline(xintercept=300, linetype="dashed", colour="grey40", linewidth=0.4) +
  labs(title="(B) Splice junction proximity — negative result",
       x="distance to nearest splice junction (bp, log)", y="density") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=9), legend.position="top")

fig5 <- (p5a | p5b) +
  plot_annotation(title="Figure 5 — Mechanism exclusion analyses",
                  theme=theme(plot.title=element_text(face="bold", size=11)))
ggsave(file.path(OUT, "Fig5_negative_results.pdf"), fig5, width=14, height=6)
message("saved Fig5")

message("\ndone. thesis figures in: ", OUT)
