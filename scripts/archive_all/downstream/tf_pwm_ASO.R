.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(PWMEnrich); library(PWMEnrich.Hsapiens.background)
  library(BSgenome.Hsapiens.UCSC.hg38); library(GenomicRanges); library(ggplot2)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/tf_motif"
KEEP_CHRS <- paste0("chr", c(1:22,"X"))
BLACKLIST <- "CENPB|ZNF274|ZNF93|UW[.]Motif"
message("loading background...")
data(PWMLogn.hg19.MotifDb.Hsap); bg <- PWMLogn.hg19.MotifDb.Hsap
gr <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
gr <- gr[gr$context=="CG" & as.character(seqnames(gr)) %in% KEEP_CHRS]
if (length(gr) > 5000) { set.seed(42); gr <- gr[sample(length(gr),5000)] }
message("n DMRs: ", length(gr))
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
res  <- motifEnrichment(seqs, bg, verbose=FALSE)
rep  <- groupReport(res)
df   <- as.data.frame(rep)
df   <- df[!grepl(BLACKLIST, df$target, ignore.case=TRUE, perl=TRUE),]
df   <- df[!duplicated(df$target),]
write.csv(df, file.path(OUT_DIR, "pwmenrich_top_motifs_ASO_5k.csv"), row.names=FALSE)
df   <- df[order(df$p.value),]
# top 20 plot
top20 <- head(df,20)
top20$target <- factor(top20$target, levels=rev(top20$target))
top20$neg_logp <- -log10(top20$p.value + 1e-300)
p <- ggplot(top20, aes(x=target, y=neg_logp)) +
  geom_col(fill="#1F3A5F") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="#C0392B", linewidth=0.5) +
  coord_flip() + theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold")) +
  labs(title=paste0("Top 20 TF motifs: ASO (n=5000)"),
       subtitle="PWMEnrich MotifDb | dashed = p=0.05",
       x=NULL, y=expression(-log[10](p)))
ggsave(file.path(OUT_DIR, "pwmenrich_report_ASO_5k.pdf"), p, width=10, height=8, device=cairo_pdf)
# top 10 plot
top10 <- head(df,10)
top10$target <- factor(top10$target, levels=rev(top10$target))
top10$neg_logp <- -log10(top10$p.value + 1e-300)
p10 <- ggplot(top10, aes(x=target, y=neg_logp)) +
  geom_col(fill="#1F3A5F") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="#C0392B", linewidth=0.5) +
  coord_flip() + theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold")) +
  labs(title=paste0("Top 10 TF motifs: ASO (n=5000)"),
       subtitle="PWMEnrich MotifDb background", x=NULL, y=expression(-log[10](p)))
ggsave(file.path(OUT_DIR, "pwmenrich_top10_ASO_5k.pdf"), p10, width=9, height=6, device=cairo_pdf)
message("Done: ASO")
