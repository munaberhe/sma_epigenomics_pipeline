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
gr <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
gr <- gr[gr$context=="CG" & as.character(seqnames(gr)) %in% KEEP_CHRS]
if (length(gr) > 5000) { set.seed(42); gr <- gr[sample(length(gr),5000)] }
message("n DMRs: ", length(gr))
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
res  <- motifEnrichment(seqs, bg, verbose=FALSE)
rep  <- groupReport(res)
# save report RDS for native plotting
saveRDS(rep, file.path(OUT_DIR, "pwmenrich_report_ASO_VPA_5k.rds"))
# filter blacklist
rep_named <- rep[!grepl(BLACKLIST, as.data.frame(rep)$target, ignore.case=TRUE, perl=TRUE),]
rep_named <- rep_named[!duplicated(as.data.frame(rep_named)$target),]
# save CSV
df <- as.data.frame(rep_named)
write.csv(df, file.path(OUT_DIR, "pwmenrich_top_motifs_ASO_VPA_5k.csv"), row.names=FALSE)
# native PWMEnrich plot - top 20
N20 <- min(20, length(rep_named))
cairo_pdf(file.path(OUT_DIR, "pwmenrich_report_ASO_VPA_5k.pdf"), width=12, height=max(10,N20*0.65))
plot(rep_named[1:N20])
dev.off()
message("saved: pwmenrich_report_ASO_VPA_5k.pdf")
# native PWMEnrich plot - top 10
N10 <- min(10, length(rep_named))
cairo_pdf(file.path(OUT_DIR, "pwmenrich_top10_ASO_VPA_5k.pdf"), width=11, height=max(8,N10*0.70))
plot(rep_named[1:N10])
dev.off()
message("Done: ASO_VPA")
