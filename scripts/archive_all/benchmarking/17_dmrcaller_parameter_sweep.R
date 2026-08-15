.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork); library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/dmrcaller_roc"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",     rds="results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds"),
  list(name="Scramble_VPA_vs_Scramble_CTRL", rds="results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds"),
  list(name="ASO_VPA_vs_Scramble_CTRL",      rds="results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
)
THRESHOLDS <- seq(0.05, 0.50, by=0.05)

sweep_one <- function(name, rds) {
  if (!file.exists(rds)) return(NULL)
  gr <- readRDS(rds)
  gr <- gr[gr$context=="CG"]
  gr$abs_diff <- abs(gr$proportion1 - gr$proportion2)
  data.frame(
    contrast=name,
    threshold=THRESHOLDS,
    n_hyper=sapply(THRESHOLDS, function(t) sum(gr$abs_diff>=t & gr$regionType=="loss")),
    n_hypo =sapply(THRESHOLDS, function(t) sum(gr$abs_diff>=t & gr$regionType=="gain")),
    total  =sapply(THRESHOLDS, function(t) sum(gr$abs_diff>=t))
  )
}

all_rows <- list()
for (ct in CONTRASTS) {
  r <- sweep_one(ct$name, ct$rds)
  if (!is.null(r)) all_rows[[ct$name]] <- r
}
big_df <- do.call(rbind, all_rows)
write.csv(big_df, file.path(OUT_DIR,"parameter_sweep_table.csv"), row.names=FALSE)

p <- ggplot(big_df, aes(x=threshold, y=total+1, colour=contrast)) +
  geom_line(linewidth=1) + geom_point(size=2) +
  scale_y_log10(labels=scales::comma) +
  scale_x_continuous(labels=scales::percent_format(1)) +
  scale_colour_manual(values=c("#1F3A5F","#C0392B","#D4A017"), name=NULL) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="top") +
  labs(title="DMR count vs methylation-difference threshold",
       subtitle="Sensitivity of DMR calls to proportion-difference cutoff",
       x="Absolute methylation difference threshold",
       y="DMRs retained (log scale)")
ggsave(file.path(OUT_DIR,"dmr_threshold_sweep.pdf"),
       p, width=8, height=5, device=cairo_pdf)
message("Done. Outputs in: ", OUT_DIR)
