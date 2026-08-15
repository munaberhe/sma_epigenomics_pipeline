.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork); library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/manhattan"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",      rds="results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",  rds="results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds"),
  list(name="ASO_VPA_vs_Scramble_CTRL",       rds="results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds"),
  list(name="ASO_VPA_vs_ASO_CTRL",            rds="results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds")
)

CHR_ORDER <- paste0("chr", c(1:22,"X"))
HG38_LEN <- c(248956422,242193529,198295559,190214555,181538259,170805979,
               159345973,145138636,138394717,133797422,135086622,133275309,
               114364328,107043718,101991189, 90338345, 83257441, 80373285,
                58617616, 64444167, 46709983, 50818468,156040895)
names(HG38_LEN) <- CHR_ORDER
offsets <- cumsum(c(0, HG38_LEN[-length(HG38_LEN)])); names(offsets) <- CHR_ORDER
mids <- offsets + HG38_LEN/2

render_manhattan <- function(name, rds) {
  if (!file.exists(rds)) { message("missing: ", rds); return(NULL) }
  gr <- readRDS(rds)
  gr <- gr[as.character(seqnames(gr)) %in% CHR_ORDER & gr$context=="CG"]
  d <- data.frame(chr=as.character(seqnames(gr)), start=start(gr),
                  p=gr$pValue, type=gr$regionType)
  d <- d[!is.na(d$p) & d$p>0,]
  d$chr <- factor(d$chr, levels=CHR_ORDER)
  d$bp_cum <- offsets[as.character(d$chr)] + d$start
  d$neg_logp <- pmin(-log10(d$p), 50)
  thr <- -log10(0.05/nrow(d))
  ggplot(d, aes(x=bp_cum, y=neg_logp, colour=chr)) +
    geom_point(alpha=0.5, size=0.4) +
    geom_hline(yintercept=thr, linetype="dashed", colour="#C0392B", linewidth=0.5) +
    scale_colour_manual(values=rep(c("#1F3A5F","#6B7280"), length.out=length(CHR_ORDER)),
                        guide="none") +
    scale_x_continuous(breaks=mids, labels=sub("chr","",names(mids)), expand=c(0.01,0)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,NA)) +
    theme_classic(base_size=11) +
    theme(panel.grid.major.y=element_line(colour="grey92"),
          plot.title=element_text(face="bold",size=12),
          axis.text.x=element_text(size=7)) +
    labs(title=paste0("Manhattan: ", gsub("_"," ",name)),
         x="Chromosome", y=expression(-log[10](p)))
}

plots <- list()
for (ct in CONTRASTS) {
  p <- render_manhattan(ct$name, ct$rds)
  if (!is.null(p)) {
    ggsave(file.path(OUT_DIR, paste0("manhattan_",ct$name,".pdf")),
           p, width=12, height=4, device=cairo_pdf)
    message("saved: manhattan_", ct$name)
    plots[[ct$name]] <- p
  }
}

if (length(plots) > 1) {
  combined <- wrap_plots(plots, ncol=1) +
    plot_annotation(title="Genome-wide DMR significance",
                    theme=theme(plot.title=element_text(face="bold",size=13)))
  ggsave(file.path(OUT_DIR,"manhattan_combined.pdf"),
         combined, width=12, height=4*length(plots), device=cairo_pdf)
  message("saved: manhattan_combined.pdf")
}
message("Done. Outputs in: ", OUT_DIR)
