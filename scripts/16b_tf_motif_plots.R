#!/usr/bin/env Rscript
# 16b_tf_motif_plots.R
# Extract and plot results from saved monaLisa motif enrichment RDS

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggplot2)
})
.libPaths(c("~/R/library", .libPaths()))

OUT_DIR <- "results/tf_motif"

message("Loading motif enrichment results...")
se <- readRDS(file.path(OUT_DIR, "motif_enrichment_results.rds"))
message("  ", nrow(se), " motifs x ", ncol(se), " bins")

# Extract results -- use as.numeric() and as.character() to avoid type issues
log_odds    <- as.numeric(assay(se, "log2enr")[, "ASO_DMR"])
neg_log_p   <- as.numeric(assay(se, "negLog10Padj")[, "ASO_DMR"])
padj        <- 10^(-neg_log_p)
motif_names <- as.character(rowData(se)[["motif.name"]])

results_df <- data.frame(
  motif_id   = rownames(se),
  motif_name = motif_names,
  log2_enr   = log_odds,
  negLog10padj = neg_log_p,
  padj       = padj,
  stringsAsFactors = FALSE
)

# Sort by enrichment
results_df <- results_df[order(results_df$log2_enr, decreasing=TRUE), ]

message("\nTop 20 enriched motifs (nominal, no significance threshold):")
print(head(results_df[, c("motif_name","log2_enr","padj")], 20), row.names=FALSE)

message("\nSummary:")
message("  Minimum padj: ", round(min(results_df$padj, na.rm=TRUE), 4))
message("  Motifs padj < 0.05: ", sum(results_df$padj < 0.05, na.rm=TRUE))
message("  Motifs padj < 0.1:  ", sum(results_df$padj < 0.1, na.rm=TRUE))
message("  Motifs padj < 0.5:  ", sum(results_df$padj < 0.5, na.rm=TRUE))

# Save full results table
write.csv(results_df,
          file.path(OUT_DIR, "motif_enrichment_all.csv"),
          row.names=FALSE)
message("Saved: motif_enrichment_all.csv")

# Check neural TFs specifically
neural_tfs <- c("ISL1","ISL2","ROBO","NEUROD1","NEUROD2","ASCL1",
                "NHLH1","NHLH2","OLIG2","NKX2","NKX6","MNX1",
                "HB9","CHAT","SOX10","PAX6","ATOH1")

message("\nNeural TF motifs:")
for (tf in neural_tfs) {
  hits <- results_df[grepl(tf, results_df$motif_name, ignore.case=TRUE), ]
  if (nrow(hits) > 0) {
    cat(sprintf("  %-12s log2enr=%+.3f  padj=%.3e\n",
                tf, hits$log2_enr[1], hits$padj[1]))
  }
}

# Plot top 20 nominal hits (no significance filter -- honest about negative result)
top_df <- head(results_df, 20)
top_df$motif_name <- factor(top_df$motif_name, levels=rev(top_df$motif_name))

p <- ggplot(top_df, aes(x=log2_enr, y=motif_name)) +
  geom_col(fill="#4b9eff", alpha=0.8) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey30") +
  labs(
    title="TF motif enrichment at ASO-specific DMRs (top 20 nominal hits)",
    subtitle=paste0(
      "151 ASO-specific DMRs vs 5x matched background\n",
      "JASPAR2020 vertebrate motifs (n=746) — no motifs reach padj<0.05\n",
      "Minimum padj = 0.18 (NFYA). Neural TFs ISL1/NeuroD1 not enriched."),
    x="Log2 enrichment (ASO DMRs vs background)",
    y=NULL,
    caption="Negative result: ASO-driven methylation changes are not at canonical TF binding sites"
  ) +
  theme_classic(base_size=11) +
  theme(
    plot.title    = element_text(face="bold"),
    plot.subtitle = element_text(size=8, colour="grey40"),
    plot.caption  = element_text(size=8, colour="#D94F3D", face="italic"),
    axis.text.y   = element_text(size=9)
  )

ggsave(file.path(OUT_DIR, "motif_enrichment_top20_nominal.pdf"),
       p, width=10, height=7)
message("Saved: motif_enrichment_top20_nominal.pdf")

# Volcano plot
results_df$highlight <- ifelse(
  grepl("ISL1|ISL2|NEUROD|ROBO|ASCL|OLIG|NKX|MNX|SOX10",
        results_df$motif_name, ignore.case=TRUE),
  "neural", "other")

p2 <- ggplot(results_df, aes(x=log2_enr, y=negLog10padj)) +
  geom_point(aes(colour=highlight), size=1.5, alpha=0.7) +
  scale_colour_manual(
    values=c("neural"="#2E9B6F", "other"="#cccccc"),
    name=NULL,
    labels=c("neural"="Neural TFs", "other"="Other")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed",
             colour="grey50", linewidth=0.5) +
  annotate("text", x=max(results_df$log2_enr)*0.8,
           y=-log10(0.05)+0.05, label="padj=0.05",
           size=3, colour="grey50") +
  geom_text(
    data=head(results_df, 5),
    aes(label=motif_name), size=2.8, hjust=-0.1, colour="#1D6FA4") +
  labs(
    title="TF motif enrichment — volcano plot",
    subtitle="746 JASPAR2020 vertebrate motifs. No motifs reach padj<0.05.",
    x="Log2 enrichment (ASO DMRs vs background)",
    y="-log10 adjusted p-value"
  ) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"))

ggsave(file.path(OUT_DIR, "motif_enrichment_volcano.pdf"),
       p2, width=10, height=7)
message("Saved: motif_enrichment_volcano.pdf")

message("\nDone. Outputs in: ", OUT_DIR)
