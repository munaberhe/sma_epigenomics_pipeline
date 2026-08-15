.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(PWMEnrich)
  library(PWMEnrich.Hsapiens.background)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/tf_motif"

# Load background once
message("loading PWMEnrich background...")
data(PWMLogn.hg19.MotifDb.Hsap)
bg <- PWMLogn.hg19.MotifDb.Hsap
KEEP_CHRS <- paste0("chr", c(1:22,"X"))
BLACKLIST  <- c("CENPB","ZNF274","ZNF93")

# Load the 3 main contrasts
gr_aso <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
gr_vpa <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
gr_combo <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
for (gr in list(gr_aso, gr_vpa, gr_combo))
  gr <- gr[gr$context=="CG" & as.character(seqnames(gr)) %in% KEEP_CHRS]

# Synergy-only: in combo but NOT in ASO alone AND NOT in VPA alone
gr_combo_f <- gr_combo[gr_combo$context=="CG" &
              as.character(seqnames(gr_combo)) %in% KEEP_CHRS]
gr_aso_f   <- gr_aso[gr_aso$context=="CG"]
gr_vpa_f   <- gr_vpa[gr_vpa$context=="CG"]

in_aso <- countOverlaps(gr_combo_f, gr_aso_f) > 0
in_vpa <- countOverlaps(gr_combo_f, gr_vpa_f) > 0
gr_synergy <- gr_combo_f[!in_aso & !in_vpa]
message("Synergy-only DMRs: ", length(gr_synergy))

# chr17 CACNG cluster DMRs
gr_cacng <- gr_combo_f[as.character(seqnames(gr_combo_f))=="chr17" &
                        start(gr_combo_f)>=66800000 &
                        end(gr_combo_f)<=67100000]
message("chr17 CACNG DMRs: ", length(gr_cacng))

SETS <- list(
  synergy_only = gr_synergy,
  chr17_cacng  = gr_cacng
)

run_pwmenrich <- function(name, gr) {
  message("\n=== ", name, " ===")
  if (length(gr) < 10) { message("  too few DMRs, skipping"); return(NULL) }
  if (length(gr) > 5000) { set.seed(42); gr <- gr[sample(length(gr),5000)] }
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  res  <- motifEnrichment(seqs, bg, verbose=FALSE)
  rep  <- groupReport(res)
  df   <- as.data.frame(rep)
  # remove blacklist
  df <- df[!sapply(df$target, function(t)
    any(sapply(BLACKLIST, function(b) grepl(b, t, ignore.case=TRUE)))),]
  write.csv(df, file.path(OUT_DIR, paste0("pwmenrich_top_motifs_",name,".csv")),
            row.names=FALSE)
  message("  saved: ", name, " (", nrow(df), " motifs, top: ",
          head(df$target[order(df$p.value)],3), ")")
}

for (nm in names(SETS)) run_pwmenrich(nm, SETS[[nm]])
message("\nDone. Outputs in: ", OUT_DIR)
