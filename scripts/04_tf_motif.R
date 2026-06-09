#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(PWMEnrich)
  library(PWMEnrich.Hsapiens.background)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(ggplot2)
  library(ggseqlogo)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/tf_motif"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

message("loading ASO-specific DMRs...")
dmr_aso <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
dmr_aso <- dmr_aso[dmr_aso$context == "CG"]
dmr_aso <- dmr_aso[as.character(seqnames(dmr_aso)) %in%
                   paste0("chr", c(1:22,"X"))]
message("  n = ", length(dmr_aso))

message("extracting sequences...")
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, dmr_aso)
message("  sequences extracted: ", length(seqs))

message("loading PWM background...")
data(PWMLogn.hg19.MotifDb.Hsap)

message("running motif enrichment (this takes a while)...")
res <- motifEnrichment(seqs, PWMLogn.hg19.MotifDb.Hsap)
saveRDS(res, file.path(OUT_DIR, "pwmenrich_results.rds"))
message("saved: pwmenrich_results.rds")

report <- groupReport(res)
saveRDS(report, file.path(OUT_DIR, "pwmenrich_report.rds"))

report_df <- as.data.frame(report)
write.csv(report_df, file.path(OUT_DIR, "pwmenrich_top_motifs.csv"),
          row.names=FALSE)
message("Top 20 motifs:")
print(head(report_df[, c("target","raw.score","p.value")], 20),
      row.names=FALSE)

# Plot top 10 as sequenceReport style
pdf(file.path(OUT_DIR, "pwmenrich_top10_report.pdf"),
    width=12, height=8)
plot(report[1:10])
dev.off()
message("saved: pwmenrich_top10_report.pdf")

# Volcano-style plot
report_df$neg_log_p <- -log10(report_df$p.value + 1e-10)
neural_tfs <- c("ISL1","ISL2","NEUROD","ROBO","ASCL","OLIG","NKX","MNX",
                "SOX10","PAX6","NHLH","HB9","CHAT")
report_df$highlight <- sapply(report_df$target, function(t)
  any(sapply(neural_tfs, function(tf) grepl(tf, t, ignore.case=TRUE))))

p_vol <- ggplot(report_df, aes(x=raw.score, y=neg_log_p)) +
  geom_point(aes(colour=highlight), size=1.5, alpha=0.7) +
  scale_colour_manual(values=c("FALSE"="#cccccc","TRUE"="#2E9B6F"),
                      labels=c("other","neural TFs")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed",
             colour="grey50", linewidth=0.5) +
  geom_text(data=head(report_df[order(report_df$p.value),], 5),
            aes(label=target), size=2.8, hjust=-0.1, colour="#1D6FA4") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold")) +
  labs(x="Raw enrichment score", y="-log10(p-value)",
       colour=NULL)

ggsave(file.path(OUT_DIR, "pwmenrich_volcano.pdf"),
       p_vol, width=10, height=7, device=cairo_pdf)
message("saved: pwmenrich_volcano.pdf")
message("done. outputs in: ", OUT_DIR)
