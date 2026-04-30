.libPaths("~/R/library")
# cache_data.R — load chr1 BSseq objects once, save as RDS for reuse.
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)

COV_DIR  <- "results/alignments/bs/by_chr"
OUT_DIR  <- "results/dmr_benchmark"
CHROM    <- "chr1"
RDS_PATH <- file.path(OUT_DIR, "chr1_data.rds")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (file.exists(RDS_PATH)) {
  message("Cache already exists at ", RDS_PATH)
  message("Delete it to force a refresh.")
  quit(save = "no", status = 0)
}

load_group_chr1 <- function(samples) {
  paths <- file.path(COV_DIR, paste0(samples, "_", CHROM, ".CpG_report.txt.gz"))
  paths <- paths[file.exists(paths)]
  message("  Loading ", length(paths), " files...")
  readBismarkPool(paths)
}

t0 <- Sys.time()
message("Loading chr1 data...")
aso_vpa  <- load_group_chr1(c("ASO_VPA_1",  "ASO_VPA_2",  "ASO_VPA_3"))
aso_ctrl <- load_group_chr1(c("ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3"))
message("  ASO_VPA CpGs:  ", length(aso_vpa))
message("  ASO_CTRL CpGs: ", length(aso_ctrl))
message("Load took: ", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), " min")

message("Saving to ", RDS_PATH, "...")
saveRDS(list(aso_vpa = aso_vpa, aso_ctrl = aso_ctrl), RDS_PATH, compress = "xz")
message("Done. Size: ", round(file.info(RDS_PATH)$size / 1e9, 2), " GB")
