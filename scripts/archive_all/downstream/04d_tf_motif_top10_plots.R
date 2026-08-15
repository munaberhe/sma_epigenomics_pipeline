.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(PWMEnrich)
  library(PWMEnrich.Hsapiens.background)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/tf_motif"
KEEP_CHRS <- paste0("chr", c(1:22,"X"))
BLACKLIST  <- "CENPB|ZNF274|ZNF93"

message("loading background...")
data(PWMLogn.hg19.MotifDb.Hsap)
bg <- PWMLogn.hg19.MotifDb.Hsap

CONTRASTS <- list(

  list(name="VPA",          rds="results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds"),
  list(name="ASO_VPA",      rds="results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds"),
  list(name="synergy_only", rds=NULL)  # loaded from pre-computed CSV
)

# For synergy_only build GRanges from the shared dip CSV
gr_combo <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
gr_combo <- gr_combo[gr_combo$context=="CG" &
              as.character(seqnames(gr_combo)) %in% KEEP_CHRS]
gr_aso   <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
gr_aso   <- gr_aso[gr_aso$context=="CG"]
gr_vpa   <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
gr_vpa   <- gr_vpa[gr_vpa$context=="CG"]
gr_synergy <- gr_combo[countOverlaps(gr_combo,gr_aso)==0 &
                       countOverlaps(gr_combo,gr_vpa)==0]
message("synergy DMRs: ", length(gr_synergy))

for (ct in CONTRASTS) {
  message("\n=== ", ct$name, " ===")

  if (ct$name == "synergy_only") {
    gr <- gr_synergy
  } else {
    gr <- readRDS(ct$rds)
    gr <- gr[gr$context=="CG" & as.character(seqnames(gr)) %in% KEEP_CHRS]
  }

  if (length(gr) > 5000) { set.seed(42); gr <- gr[sample(length(gr),5000)] }
  message("  n DMRs: ", length(gr))

  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  res  <- motifEnrichment(seqs, bg, verbose=FALSE)
  rep  <- groupReport(res)

  # save PWMEnrich native plot (top 20)
  fname <- file.path(OUT_DIR, paste0("pwmenrich_report_",ct$name,".pdf"))
  cairo_pdf(fname, width=12, height=8)
  plot(rep[1:20])
  dev.off()
  message("saved: ", basename(fname))

  # also save top 10 as cleaner PDF
  fname2 <- file.path(OUT_DIR, paste0("pwmenrich_top10_",ct$name,".pdf"))
  cairo_pdf(fname2, width=10, height=6)
  plot(rep[1:10])
  dev.off()
  message("saved: ", basename(fname2))
}
message("\nDone.")
