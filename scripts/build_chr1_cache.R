.libPaths("~/R/library")
library(DMRcaller)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/dmr_benchmark"
CHROM   <- "chr1"
RDS_PATH <- file.path(OUT_DIR, "chr1_data.rds")

message("Building chr1 cache...")
load_one <- function(samples) {
  paths <- file.path(COV_DIR, paste0(samples, "_", CHROM, ".CpG_report.txt.gz"))
  paths <- paths[file.exists(paths)]
  message("  pooling: ", paste(basename(paths), collapse=", "))
  readBismarkPool(paths)
}

aso_vpa  <- load_one(paste0("ASO_VPA_",  1:3))
aso_ctrl <- load_one(paste0("ASO_CTRL_", 1:3))

message("ASO_VPA CpGs:  ", format(length(aso_vpa),  big.mark=","))
message("ASO_CTRL CpGs: ", format(length(aso_ctrl), big.mark=","))

saveRDS(list(aso_vpa=aso_vpa, aso_ctrl=aso_ctrl), RDS_PATH)
message("Saved cache: ", RDS_PATH)
