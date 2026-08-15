#!/usr/bin/env Rscript
.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

COLS <- c(
  ASO_CTRL      = "#1B4F8A",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#6B7280"
)

setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/tss_metaplot'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)


# DMR heatmap
message("Loading cache and DMR results...")
cache   <- readRDS('results/dmr/meth_pooled_cache.rds')
# Filter cache to CG context only
cache <- lapply(cache, function(m) m[mcols(m)$context=="CG"])
dmr_aso <- readRDS('results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds')
# Filter DMRs to CG context
dmr_aso <- dmr_aso[mcols(dmr_aso)$context=="CG"]
message("  DMRs (CG, ASO_VPA vs Scramble_CTRL): ", length(dmr_aso))
message("Building DMR heatmap...")
meth_diff <- abs(mcols(dmr_aso)$proportion1 - mcols(dmr_aso)$proportion2)
dmr_top <- dmr_aso[order(meth_diff, decreasing=TRUE)[1:500]]
message("  Top 500 by meth diff, min diff: ", round(min(sort(meth_diff, decreasing=TRUE)[1:500]), 3))

get_meth_at_dmrs <- function(meth, dmrs) {
  ov <- findOverlaps(dmrs, meth, ignore.strand=TRUE)
  vals <- rep(NA, length(dmrs))
  for (i in seq_along(dmrs)) {
    hits <- subjectHits(ov[queryHits(ov)==i])
    if (length(hits) > 0) {
      m <- sum(mcols(meth)$readsM[hits])
      n <- sum(mcols(meth)$readsN[hits])
      if (n > 0) vals[i] <- m/n
    }
  }
  vals
}

heatmap_mat <- sapply(names(cache), function(cond) {
  message("  extracting methylation for ", cond)
  get_meth_at_dmrs(cache[[cond]], dmr_top)
})
colnames(heatmap_mat) <- names(cache)

# Remove rows with too many NAs
keep <- rowSums(is.na(heatmap_mat)) < 2
heatmap_mat <- heatmap_mat[keep, ]
message("  DMRs with data: ", nrow(heatmap_mat))

col_fun <- colorRamp2(c(0, 0.5, 1),
                      c("#2166AC", "white", "#B2182B"))
col_annot <- HeatmapAnnotation(
  condition = colnames(heatmap_mat),
  col = list(condition = COLS),
  annotation_name_side = "left"
)

pdf(file.path(OUT, "DMR_heatmap_top500_ASO_VPA_methdiff.pdf"), width=8, height=12)
Heatmap(heatmap_mat,
        name="CpG meth",
        col=col_fun,
        top_annotation=col_annot,
        show_row_names=FALSE,
        show_column_names=TRUE,
        cluster_rows=TRUE,
        cluster_columns=TRUE,
        column_title="Top 500 ASO DMRs — sample clustering",
        row_title=paste0(nrow(heatmap_mat), " DMRs"),
        use_raster=TRUE)
dev.off()
message("saved: DMR_heatmap_top500.pdf")

message("done. outputs in: ", OUT)
