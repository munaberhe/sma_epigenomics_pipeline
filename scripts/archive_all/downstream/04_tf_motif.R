#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(PWMEnrich)
  library(PWMEnrich.Hsapiens.background)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(ggplot2)
  library(ggrepel)
})
source("scripts/pipeline/00_sma_palette.R")
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

# ---- TF motif enrichment at ASO DMR sequences ----
# PWMEnrich against MotifDb human PWMs.
# Caveat: PWMLogn.hg19.MotifDb.Hsap background is hg19; DMRs are hg38. Should
# match assemblies eventually, but result is "no enrichment" either way.

OUT_DIR <- "results/tf_motif"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ---- input DMRs ----
message("loading ASO DMRs...")
dmr_aso <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
dmr_aso <- dmr_aso[dmr_aso$context == "CG"]
dmr_aso <- dmr_aso[as.character(seqnames(dmr_aso)) %in%
                   paste0("chr", c(1:22,"X"))]
message("  n DMRs: ", length(dmr_aso))

message("extracting sequences...")
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, dmr_aso)
message("  sequences extracted: ", length(seqs))

# ---- enrichment ----
message("loading PWM background...")
data(PWMLogn.hg19.MotifDb.Hsap)

message("running motif enrichment (slow)...")
res <- motifEnrichment(seqs, PWMLogn.hg19.MotifDb.Hsap)
saveRDS(res, file.path(OUT_DIR, "pwmenrich_results.rds"))

report <- groupReport(res)
saveRDS(report, file.path(OUT_DIR, "pwmenrich_report.rds"))

report_df <- as.data.frame(report)
# filter out unannotated MotifDb entries (UW.Motif.xxxx have no TF name)
write.csv(report_df, file.path(OUT_DIR, "pwmenrich_top_motifs.csv"),
          row.names=FALSE)
message("top 20 motifs:")
print(head(report_df[, c("target","raw.score","p.value")], 20),
      row.names=FALSE)

# use named-TF-only report (filter UW.Motif.xxxx unannotated entries)
report_named <- report[!grepl("^UW\\.Motif\\.", as.data.frame(report)$target), ]
N <- min(30, length(report_named))
pdf(file.path(OUT_DIR, "pwmenrich_top30_report.pdf"), width=11, height=max(8, N * 0.35))
plot(report_named[1:N])

# ---- volcano ----
# nothing crosses p=0.05 in current run; volcano still useful to show
# distribution + top hits. Use ggrepel to keep labels readable.
report_df$neg_log_p <- -log10(pmax(report_df$p.value, 1e-300))
report_df$sig_level <- ifelse(report_df$p.value < 0.01, "top",
                       ifelse(report_df$p.value < 0.05, "significant", "ns"))
report_df <- report_df[order(report_df$sig_level == "ns", decreasing=TRUE), ]
report_df_named <- report_df[!grepl("^UW[.]Motif[.]", report_df$target) & !grepl("CENPB", report_df$target, ignore.case=TRUE), ]
label_df <- report_df_named[report_df_named$p.value < 0.05, ]

# soft-flag SMA-relevant TFs: motor neuron + general neural development


p_vol <- ggplot(report_df, aes(x=raw.score, y=neg_log_p)) +
  geom_point(aes(colour=sig_level), size=1.6, alpha=0.7) +
  scale_colour_manual(values=c("ns"="#cccccc","significant"="#2E9B6F","top"="#C0392B"),
                     labels=c("ns","p<0.05","p<0.01"), name=NULL) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed",
             colour="grey50", linewidth=0.5) +
  geom_text_repel(data=label_df, aes(x=raw.score, label=target, colour=sig_level),
                  size=3, max.overlaps=Inf, min.segment.length=0,
                  segment.colour="grey60", segment.size=0.3,
                  box.padding=0.4, point.padding=0.3,
                  show.legend=FALSE) +
  annotate("text", x=Inf, y=-log10(0.05), hjust=1.05, vjust=-0.5,
           label="p = 0.05", size=2.8, colour="grey40") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"),
        legend.position="top") +
  scale_x_log10() +
  labs(title="TF motif enrichment at ASO DMRs (top hits labeled)",
       subtitle=paste0("PWMEnrich MotifDb.Hsap (x-axis log10), n=", length(dmr_aso),
                       " DMRs | green=p<0.05, red=p<0.01 | all significant TFs labelled"),
       x="Raw enrichment score",
       y="-log10(p-value)")

ggsave(file.path(OUT_DIR, "pwmenrich_volcano.pdf"),
       p_vol, width=10, height=7, device=cairo_pdf)
message("saved: pwmenrich_volcano.pdf")

# ---- top-N table (slide-ready, since nothing passes 0.05) ----
top_n <- 15
top_table <- head(report_df[order(report_df$p.value),
                            c("target","raw.score","p.value")], top_n)
top_table$p.value  <- signif(top_table$p.value, 3)
top_table$raw.score <- signif(top_table$raw.score, 3)
write.table(top_table,
            file.path(OUT_DIR, "pwmenrich_top15_table.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
message("saved: pwmenrich_top15_table.tsv")

message("done. outputs in: ", OUT_DIR)
