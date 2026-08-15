.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/h3k9me2_overlap"

# ---- reload the already-computed outputs from 03a_h3k9me2.R ----
plot_df <- readRDS(file.path(OUT_DIR, "h3k9me2_plot_df.rds"))
summary_tsv <- read.table(file.path(OUT_DIR, "h3k9me2_signal_summary.tsv"),
                          sep="\t", header=TRUE)
smn2_df <- read.table(file.path(OUT_DIR, "smn2_h3k9me2_signal.tsv"),
                      sep="\t", header=TRUE)

# ---- recompute delta values from plot_df source data is not directly
# available here, so re-derive what we need from the saved summary and
# from a fresh read of the underlying per-region signal if present.
# Simpler: reload sig_aso / sig_bg_aso directly if cached, else recompute
# minimal delta summary from h3k9me2_signal_summary.tsv (CTR-only) is not
# enough for delta, so we pull delta stats from the already-saved
# h3k9me2_delta_signal.pdf run's printed values instead -- re-derive here
# by re-reading the same RDS DMR/background signal objects if cached.
# For robustness, recompute delta directly from plot_df2-equivalent source:
# plot_df (CTR vs BG only) doesn't carry ASO arm, so fall back to rerunning
# the delta computation lightly using the saved per-region tables if they
# exist; otherwise this section is skipped gracefully.

delta_cache <- file.path(OUT_DIR, "h3k9me2_delta_cache.rds")
if (file.exists(delta_cache)) {
  delta_df <- readRDS(delta_cache)
} else {
  message("No cached delta table found -- run 03a_h3k9me2.R first to generate it.")
  message("Falling back to summary-only panel B for now.")
  delta_df <- NULL
}

# ---- Panel A: genome-wide delta H3K9me2, ASO DMRs vs background ----
if (!is.null(delta_df)) {
  wtest <- wilcox.test(delta_df$delta[delta_df$group == "ASO DMRs"],
                       delta_df$delta[grepl("Background", delta_df$group)],
                       alternative = "two.sided")
  n_dmr <- sum(delta_df$group == "ASO DMRs")
  n_bg  <- sum(grepl("Background", delta_df$group))

  panel_a <- ggplot(delta_df, aes(x=group, y=delta, fill=group)) +
    geom_boxplot(outlier.shape=NA, width=0.4) +
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="white", colour="#1A1A1A", stroke=0.8) +
    stat_summary(fun=mean, geom="text", aes(label=sprintf("mean = %.2f", after_stat(y))), hjust=-0.18, vjust=0.4, size=4.2, fontface="bold", colour="#1A1A1A") +
    geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    scale_fill_manual(values=c("Background"="#cccccc", "ASO DMRs"="#1B4F8A"),
                      guide="none") +
    coord_cartesian(ylim=c(-15, 15)) +
    annotate("text", x=1.5, y=14,
             label=sprintf("Wilcoxon p = %.2e", wtest$p.value),
             size=3.6, colour="grey20") +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold", size=12),
          plot.subtitle=element_text(size=9.5, colour="grey30")) +
    labs(title="B. Genome-wide: ASO DMRs vs matched background",
         subtitle=sprintf("n=%d ASO-defined WGBS DMRs, n=%d chr/width-matched background regions\nMean delta: DMRs=%.2f, background=%.2f | H3K9me2 ChIP-seq, Marasco et al. 2022 (GSE167762)", n_dmr, n_bg, mean(delta_df$delta[delta_df$group=="ASO DMRs"]), mean(delta_df$delta[grepl("Background", delta_df$group)])),
         x=NULL, y="Delta H3K9me2 signal (ASO - CTR)")
} else {
  panel_a <- ggplot() + theme_void() +
    labs(title="Panel A unavailable: rerun 03a_h3k9me2.R to cache delta values")
}

# ---- Panel B: SMN2-specific signal, CTR vs ASO ----
smn2_df$arm <- ifelse(grepl("^CTR", smn2_df$condition), "CTR", "ASO")
smn2_summary <- aggregate(signal ~ arm, data=smn2_df, FUN=mean)
smn2_summary$arm <- factor(smn2_summary$arm, levels=c("CTR","ASO"))
fold <- smn2_summary$signal[smn2_summary$arm=="ASO"] /
        smn2_summary$signal[smn2_summary$arm=="CTR"]

panel_b <- ggplot(smn2_summary, aes(x=arm, y=signal, fill=arm)) +
  geom_col(width=0.5, colour="grey20") +
  geom_point(data=smn2_df, aes(x=arm, y=signal), inherit.aes=FALSE,
             size=2, shape=21, fill="white", colour="grey30") +
  scale_fill_manual(values=c("CTR"="#6B7280","ASO"="#1B4F8A"), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  geom_text(aes(label=sprintf("%.2f", signal)), data=smn2_summary, vjust=-0.6, size=5, fontface="bold", colour="#1A1A1A") +
  annotate("text", x=1.5, y=max(smn2_df$signal)*1.05,
           label=sprintf("%.2fx increase", fold), size=3.8, fontface="bold") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold", size=12),
        plot.subtitle=element_text(size=9.5, colour="grey30")) +
  labs(title="A. SMN2 locus: H3K9me2 signal by condition",
       subtitle=sprintf("CTR=%.2f, ASO=%.2f, chr5:70,049,638-70,078,522\nReplicate points shown as open circles",
                        smn2_summary$signal[smn2_summary$arm=="CTR"],
                        smn2_summary$signal[smn2_summary$arm=="ASO"]),
       x=NULL, y="H3K9me2 ChIP-seq signal")

combined <- panel_b | panel_a
combined <- combined +
  plot_annotation(
    title="H3K9me2 at the SMN2 locus is elevated under ASO treatment, but ASO-defined DMRs show no genome-wide H3K9me2 enrichment",
    theme=theme(plot.title=element_text(face="bold", size=12.5))
  )

ggsave(file.path(OUT_DIR, "h3k9me2_twopanel_final.pdf"), combined,
       width=13, height=6, device=cairo_pdf)
message("saved: h3k9me2_twopanel_final.pdf")
