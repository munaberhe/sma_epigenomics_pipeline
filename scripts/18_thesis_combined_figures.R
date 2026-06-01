#!/usr/bin/env Rscript
# 18_thesis_combined_figures.R
# Generates combined multi-panel thesis figures from saved RDS data
# All figures sized and styled for thesis submission
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(SummarizedExperiment)
})
.libPaths(c("~/R/library", .libPaths()))

setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT <- "results/thesis_figures"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# Consistent colour scheme
COLS <- c(
  ASO_CTRL      = "#2E9B6F",
  ASO_VPA       = "#D94F3D",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#1D6FA4",
  Hypo          = "#4393C3",
  Hyper         = "#D6604D"
)

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       label="ASO alone vs Scramble CTRL\n(nusinersen off-target signal)"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       label="VPA alone vs Scramble CTRL\n(HDAC inhibitor effect)"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       label="ASO+VPA vs ASO alone\n(VPA effect on ASO background)")
)

CHROMS <- paste0("chr", c(1:22, "X", "Y"))

message("Loading DMR objects...")
dmrs <- lapply(CONTRASTS, function(ct) {
  d <- readRDS(paste0("results/dmr/dmr_", ct$name, ".rds"))
  d[mcols(d)$cytosinesCount >= 6]
})
names(dmrs) <- sapply(CONTRASTS, function(ct) ct$name)

# ── FIGURE 1: METHYLATION DIFFERENCE HISTOGRAMS (3 panels) ───
message("Building Figure 1: Methylation difference histograms...")

make_hist <- function(dmr_obj, label, tag) {
  df <- data.frame(
    methDiff  = mcols(dmr_obj)$proportion1 - mcols(dmr_obj)$proportion2,
    direction = ifelse(mcols(dmr_obj)$regionType == "gain", "Hypo", "Hyper")
  )
  n_total <- nrow(df)
  n_hypo  <- sum(df$direction == "Hypo")
  n_hyper <- sum(df$direction == "Hyper")

  ggplot(df, aes(x=methDiff, fill=direction)) +
    geom_histogram(bins=60, colour="white", linewidth=0.1) +
    geom_vline(xintercept=0, linewidth=0.5, linetype="dashed", colour="grey20") +
    geom_vline(xintercept=c(-0.2, 0.2), linewidth=0.3,
               linetype="dotted", colour="grey50") +
    scale_fill_manual(values=c(Hypo=unname(COLS["Hypo"]), Hyper=unname(COLS["Hyper"])),
                      name=NULL) +
    scale_x_continuous(limits=c(-1,1), breaks=seq(-1,1,0.5),
                       labels=scales::percent_format(accuracy=1)) +
    scale_y_continuous(labels=scales::comma) +
    labs(title=paste0("(", tag, ") ", label),
         x="Methylation difference (treatment − reference)",
         y="Number of DMRs",
         caption=paste0("n=", format(n_total, big.mark=","),
                        "  Hypo=", format(n_hypo, big.mark=","),
                        "  Hyper=", format(n_hyper, big.mark=","))) +
    theme_bw(base_size=10) +
    theme(plot.title=element_text(face="bold", size=9),
          plot.caption=element_text(size=7, colour="grey40"),
          legend.position="top",
          panel.grid.minor=element_blank())
}

p1a <- make_hist(dmrs[["ASO_CTRL_vs_Scramble_CTRL"]],
                 "ASO alone vs Scramble CTRL", "A")
p1b <- make_hist(dmrs[["Scramble_VPA_vs_Scramble_CTRL"]],
                 "VPA alone vs Scramble CTRL", "B")
p1c <- make_hist(dmrs[["ASO_VPA_vs_ASO_CTRL"]],
                 "ASO+VPA vs ASO alone", "C")

fig1 <- (p1a | p1b | p1c) +
  plot_annotation(
    title="Figure 1 — Methylation difference distributions across three contrasts",
    subtitle=paste0("(A) ASO produces a bidirectional signal. ",
                    "(B) VPA produces near-exclusively hypomethylated DMRs. ",
                    "(C) VPA on ASO background reproduces the VPA-alone pattern,\n",
                    "confirming that the two drugs act through independent epigenetic mechanisms."),
    theme=theme(plot.title=element_text(face="bold", size=11),
                plot.subtitle=element_text(size=8, colour="grey30"))
  )

ggsave(file.path(OUT, "Fig1_methylation_diff_histograms.pdf"),
       fig1, width=15, height=5)
ggsave(file.path(OUT, "Fig1_methylation_diff_histograms.png"),
       fig1, width=15, height=5, dpi=300)
message("Saved: Fig1_methylation_diff_histograms")

# ── FIGURE 2: PER-CHROMOSOME BARS (2 panels) ──────────────────
message("Building Figure 2: Per-chromosome bar charts...")

make_chr_bar <- function(dmr_obj, label, tag) {
  chr_counts <- as.data.frame(table(
    chr       = as.character(seqnames(dmr_obj)),
    direction = ifelse(mcols(dmr_obj)$regionType == "gain", "Hypo", "Hyper")
  ))
  chr_counts$chr <- factor(chr_counts$chr, levels=CHROMS)
  chr_counts <- chr_counts[chr_counts$chr %in% CHROMS, ]
  chr_counts$Freq_signed <- ifelse(chr_counts$direction=="Hypo",
                                   -chr_counts$Freq, chr_counts$Freq)

  ggplot(chr_counts, aes(x=chr, y=Freq_signed, fill=direction)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept=0, linewidth=0.4, colour="grey30") +
    scale_fill_manual(values=c(Hypo=unname(COLS["Hypo"]), Hyper=unname(COLS["Hyper"])),
                      name=NULL) +
    scale_y_continuous(labels=function(x) format(abs(x), big.mark=",")) +
    labs(title=paste0("(", tag, ") ", label),
         x="Chromosome", y="Number of DMRs") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
          plot.title=element_text(face="bold", size=9),
          legend.position="top",
          panel.grid.minor=element_blank())
}

p2a <- make_chr_bar(dmrs[["ASO_CTRL_vs_Scramble_CTRL"]],
                    "ASO alone vs Scramble CTRL", "A")
p2b <- make_chr_bar(dmrs[["Scramble_VPA_vs_Scramble_CTRL"]],
                    "VPA alone vs Scramble CTRL", "B")

fig2 <- (p2a / p2b) +
  plot_annotation(
    title="Figure 2 — DMR chromosomal distribution",
    subtitle=paste0("(A) ASO shows disproportionate chrX enrichment (620 DMRs, 18% of total), ",
                    "absent in (B) VPA alone.\n",
                    "VPA DMRs scale with chromosome size, consistent with non-specific HDAC inhibition."),
    theme=theme(plot.title=element_text(face="bold", size=11),
                plot.subtitle=element_text(size=8, colour="grey30"))
  )

ggsave(file.path(OUT, "Fig2_per_chromosome_bars.pdf"),
       fig2, width=14, height=10)
ggsave(file.path(OUT, "Fig2_per_chromosome_bars.png"),
       fig2, width=14, height=10, dpi=300)
message("Saved: Fig2_per_chromosome_bars")

# ── FIGURE 3: QC TRIPTYCH ─────────────────────────────────────
message("Building Figure 3: QC triptych...")

# Load correlation matrix for heatmap
cor_mat <- read.csv("results/qc/methylation/methylation_correlation.csv",
                    row.names=1)

# Melt for ggplot
cor_melt <- reshape2::melt(as.matrix(cor_mat))
cor_melt$Var1 <- factor(cor_melt$Var1)
cor_melt$Var2 <- factor(cor_melt$Var2)

p3a <- ggplot(cor_melt, aes(Var1, Var2, fill=value)) +
  geom_tile(colour="white", linewidth=0.3) +
  geom_text(aes(label=round(value,2)), size=2.2, fontface="bold") +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                       midpoint=0.87, limits=c(0.84,1),
                       name="Pearson r") +
  labs(title="(A) Sample correlation heatmap",
       x=NULL, y=NULL) +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        axis.text.y=element_text(size=7),
        plot.title=element_text(face="bold", size=9),
        panel.grid=element_blank())

# Per-condition violin from DMR proportions
dmr_aso <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
dmr_aso <- dmr_aso[mcols(dmr_aso)$cytosinesCount >= 6]
dmr_vpa <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
dmr_vpa <- dmr_vpa[mcols(dmr_vpa)$cytosinesCount >= 6]
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
  geom_boxplot(width=0.07, fill="white",
               outlier.size=0.2, outlier.alpha=0.2) +
  scale_fill_manual(values=unname(COLS[c("ASO_CTRL","Scramble_CTRL","Scramble_VPA","ASO_VPA")]), labels=c("ASO_CTRL","Scramble_CTRL","Scramble_VPA","ASO_VPA")) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(title="(B) CpG methylation distribution at DMR loci",
       x=NULL, y="CpG methylation") +
  theme_bw(base_size=9) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=30, hjust=1, size=8),
        plot.title=element_text(face="bold", size=9))

fig3 <- (p3a | p3b) +
  plot_annotation(
    title="Figure 3 — Quality control",
    subtitle=paste0("(A) Per-replicate Pearson correlation confirms within-condition clustering. ",
                    "ASO_VPA shows lower within-group correlation (r=0.845-0.858)\n",
                    "reflecting combined epigenetic variance. ",
                    "(B) VPA-treated conditions show global hypomethylation shift."),
    theme=theme(plot.title=element_text(face="bold", size=11),
                plot.subtitle=element_text(size=8, colour="grey30"))
  )

ggsave(file.path(OUT, "Fig3_QC.pdf"),   fig3, width=14, height=6)
ggsave(file.path(OUT, "Fig3_QC.png"),   fig3, width=14, height=6, dpi=300)
message("Saved: Fig3_QC")

# ── FIGURE 4: H3K9me2 ─────────────────────────────────────────
message("Building Figure 4: H3K9me2...")

# Re-read the signal data saved by script 11
smn2_df <- read.table("results/h3k9me2_overlap/smn2_h3k9me2_signal.tsv",
                      header=TRUE, sep="\t")
summary_df <- read.table("results/h3k9me2_overlap/h3k9me2_signal_summary.tsv",
                         header=TRUE, sep="\t")

# Rebuild plot_df2 from DMR objects
dmr_aso_hc <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
dmr_aso_hc <- dmr_aso_hc[mcols(dmr_aso_hc)$cytosinesCount >= 6]

# Check what columns exist in smn2_df
message("SMN2 df columns: ", paste(colnames(smn2_df), collapse=", "))
message("Summary df: ")
print(summary_df)

# SMN2 locus bar chart
if (all(c("condition","signal") %in% colnames(smn2_df))) {
  smn2_df$condition <- factor(smn2_df$condition,
    levels=c("CTR_R1","CTR_R2","ASO_R1","ASO_R2"))
  smn2_df$treatment <- ifelse(grepl("CTR", smn2_df$condition), "CTR", "ASO")

  p4b <- ggplot(smn2_df, aes(x=condition, y=signal, fill=treatment)) +
    geom_col() +
    scale_fill_manual(values=c(CTR="#1D6FA4", ASO="#2E9B6F"),
                      name="Treatment") +
    labs(title="(B) H3K9me2 signal at SMN2 locus",
         subtitle="2.2-fold increase in ASO vs CTR (Marasco et al. 2022, GSE167762)",
         x=NULL, y="Mean H3K9me2 signal") +
    theme_bw(base_size=10) +
    theme(plot.title=element_text(face="bold", size=9),
          plot.subtitle=element_text(size=8))
} else {
  p4b <- ggplot() +
    annotate("text", x=0.5, y=0.5, label="SMN2 H3K9me2\n(see smn2_h3k9me2_signal.tsv)",
             size=4) + theme_void()
}

# Genome-wide boxplot from saved RDS
plot_df <- readRDS("results/h3k9me2_overlap/h3k9me2_plot_df.rds")
plot_df$group <- factor(plot_df$group, levels=c("Background","ASO DMRs"))

p4a <- ggplot(plot_df, aes(x=group, y=signal, fill=group)) +
  geom_boxplot(outlier.size=0.3, outlier.alpha=0.3) +
  facet_wrap(~set, scales="free_y") +
  scale_fill_manual(values=c("Background"="#cccccc", "ASO DMRs"="#2E9B6F"),
                    name=NULL) +
  scale_y_continuous(limits=c(0, quantile(plot_df$signal, 0.99, na.rm=TRUE))) +
  labs(title="(A) H3K9me2 signal at ASO DMRs vs matched background",
       subtitle="No genome-wide enrichment (Wilcoxon p=1.0) — locus-specific mechanism",
       x=NULL, y="Mean H3K9me2 signal") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=8),
        legend.position="none",
        strip.text=element_text(face="bold"))

fig4 <- (p4a | p4b) +
  plot_annotation(
    title="Figure 4 — H3K9me2 enrichment analysis",
    subtitle=paste0("(A) ASO DMRs show no genome-wide H3K9me2 enrichment vs matched background ",
                    "(Wilcoxon p=1.0), confirming locus-specific mechanism.\n",
                    "(B) H3K9me2 is 2.2-fold elevated at the SMN2 locus specifically ",
                    "(CTR=7.16, ASO=15.53), validating kinetic coupling."),
    theme=theme(plot.title=element_text(face="bold", size=11),
                plot.subtitle=element_text(size=8, colour="grey30"))
  )

ggsave(file.path(OUT, "Fig4_H3K9me2.pdf"), fig4, width=12, height=5)
ggsave(file.path(OUT, "Fig4_H3K9me2.png"), fig4, width=12, height=5, dpi=300)
message("Saved: Fig4_H3K9me2")

# ── FIGURE 5: NEGATIVE RESULTS ────────────────────────────────
message("Building Figure 5: Negative results...")

# TF motif volcano
se <- readRDS("results/tf_motif/motif_enrichment_results.rds")
motif_df <- data.frame(
  motif   = rownames(se),
  log2enr = as.numeric(SummarizedExperiment::assay(se, "log2enr")[,1]),
  neglog10padj = as.numeric(SummarizedExperiment::assay(se, "negLog10Padj")[,1])
)
motif_df <- motif_df[complete.cases(motif_df), ]
motif_df$padj <- 10^(-motif_df$neglog10padj)
motif_df$significant <- motif_df$padj < 0.05
motif_df$neglog10p <- motif_df$neglog10padj

top_nominal <- head(motif_df[order(motif_df$padj),], 5)

p5a <- ggplot(motif_df, aes(x=log2enr, y=neglog10p)) +
  geom_point(aes(colour=significant), alpha=0.5, size=1) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed",
             colour="red", linewidth=0.5) +
  scale_colour_manual(values=c("FALSE"="grey70","TRUE"="red"),
                      guide="none") +
  geom_text(data=top_nominal,
            aes(label=gsub("_.*","",motif)),
            size=2.5, vjust=-0.5, colour="grey30") +
  labs(title="(A) TF motif enrichment — negative result",
       subtitle="No significant enrichment (min p.adj=0.18, 746 motifs tested)",
       x="log2 enrichment", y="-log10(p.adj)") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=8))

# Splice junction density
sj_df <- read.csv("results/splice_junction/splice_junction_distances.csv")
sj_summary <- read.csv("results/splice_junction/splice_junction_summary.csv")

dmr_med <- sj_summary$median_dist_bp[sj_summary$group=="ASO_specific_DMRs"]
bg_med  <- sj_summary$median_dist_bp[sj_summary$group=="Background"]
wtest_p <- sj_summary$wilcox_p[sj_summary$group=="ASO_specific_DMRs"]

p5b <- ggplot(sj_df, aes(x=distance, fill=type, colour=type)) +
  geom_density(alpha=0.4, linewidth=0.7) +
  scale_x_log10(labels=scales::comma) +
  scale_fill_manual(values=c("ASO DMRs"="#2E9B6F", "Background"="grey70"), name=NULL) +
  scale_colour_manual(values=c("ASO DMRs"="#2E9B6F", "Background"="grey50"), name=NULL) +
  geom_vline(xintercept=300, linetype="dashed",
             colour="grey40", linewidth=0.4) +
  labs(title="(B) Splice junction proximity — negative result",
       subtitle=paste0("DMR median=", format(dmr_med, big.mark=","),
                       "bp vs background=", format(bg_med, big.mark=","),
                       "bp  Wilcoxon p=", wtest_p),
       x="Distance to nearest splice junction (bp, log scale)",
       y="Density") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=8),
        legend.position="top")

fig5 <- (p5a | p5b) +
  plot_annotation(
    title="Figure 5 — Mechanism exclusion analyses",
    subtitle=paste0("(A) No TF binding site disruption detected (min p.adj=0.18, JASPAR2020). ",
                    "(B) ASO-specific DMRs are not enriched near splice junctions\n",
                    "(Wilcoxon p=0.939), ruling out genome-wide kinetic coupling ",
                    "as the off-target mechanism."),
    theme=theme(plot.title=element_text(face="bold", size=11),
                plot.subtitle=element_text(size=8, colour="grey30"))
  )

ggsave(file.path(OUT, "Fig5_negative_results.pdf"), fig5, width=14, height=6)
ggsave(file.path(OUT, "Fig5_negative_results.png"), fig5, width=14, height=6, dpi=300)
message("Saved: Fig5_negative_results")

message("\nAll thesis figures saved to: ", OUT)
list.files(OUT)
