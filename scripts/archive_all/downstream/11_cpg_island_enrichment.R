.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges); library(annotatr)
  library(ggplot2); library(dplyr); library(patchwork)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/cpg_context"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

message("Building CpG annotations...")
cpg_anns <- build_annotations(genome="hg38",
  annotations=c("hg38_cpg_islands","hg38_cpg_shores","hg38_cpg_shelves","hg38_cpg_inter"))
hg38_size <- 3.1e9
bg_sizes  <- tapply(width(cpg_anns), cpg_anns$type, sum)
bg_frac   <- bg_sizes / hg38_size

CONTRASTS <- list(
  list(name="ASO_specific", rds="results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds"),
  list(name="VPA",          rds="results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds"),
  list(name="ASO_VPA",      rds="results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
)
CHRS <- paste0("chr",c(1:22,"X"))

calc_enrich <- function(name, rds) {
  if (!file.exists(rds)) return(NULL)
  gr <- readRDS(rds)
  gr <- gr[as.character(seqnames(gr)) %in% CHRS & gr$context=="CG"]
  rows <- list()
  for (dir in c("loss","gain")) {
    sub <- gr[gr$regionType==dir]
    if (length(sub)==0) next
    for (cls in names(bg_frac)) {
      ann_sub <- cpg_anns[cpg_anns$type==cls]
      obs <- length(subsetByOverlaps(sub, ann_sub)) / length(sub)
      exp <- bg_frac[[cls]]
      rows[[length(rows)+1]] <- data.frame(
        contrast=name, direction=dir, context=cls,
        observed=obs, expected=exp,
        log2fe=log2((obs+1e-6)/(exp+1e-6)))
    }
  }
  do.call(rbind, rows)
}

all_rows <- list()
for (ct in CONTRASTS) {
  r <- calc_enrich(ct$name, ct$rds)
  if (!is.null(r)) {
    all_rows[[ct$name]] <- r
    write.csv(r, file.path(OUT_DIR, paste0("cpg_context_",ct$name,".csv")), row.names=FALSE)
  }
}

combined_df <- do.call(rbind, all_rows)
combined_df$context <- factor(combined_df$context,
  levels=c("hg38_cpg_islands","hg38_cpg_shores","hg38_cpg_shelves","hg38_cpg_inter"),
  labels=c("CpG island","Shore","Shelf","Open sea"))
combined_df$direction <- factor(combined_df$direction,
  levels=c("loss","gain"), labels=c("Hyper","Hypo"))

p <- ggplot(combined_df, aes(x=context, y=log2fe, fill=direction)) +
  geom_col(position=position_dodge(0.85), width=0.8) +
  geom_hline(yintercept=0, colour="grey40") +
  scale_fill_manual(values=c("Hypo"="#1F3A5F","Hyper"="#C0392B"), name=NULL) +
  facet_wrap(~contrast, ncol=1) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"),
        axis.text.x=element_text(angle=20, hjust=1)) +
  labs(title="DMR enrichment in CpG-island context",
       x=NULL, y=expression(log[2]~"(observed/expected)"))

ggsave(file.path(OUT_DIR,"cpg_context_enrichment.pdf"),
       p, width=7, height=3.5*length(all_rows), device=cairo_pdf)
message("Done. Outputs in: ", OUT_DIR)
