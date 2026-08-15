.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges); library(ggplot2); library(dplyr)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/meth_distribution"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

CONTRASTS <- list(
  list(name="ASO_CTRL",  rds="results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds",  cond="ASO_CTRL"),
  list(name="VPA",       rds="results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds", cond="Scramble_VPA"),
  list(name="ASO_VPA",   rds="results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds",   cond="ASO_VPA")
)

load_cx_chr5 <- function(condition) {
  grs <- lapply(1:3, function(r) {
    path <- file.path("results/alignments_smn1_masked/chr5_cx",
                      sprintf("%s_%d_chr5.CX_report.txt", condition, r))
    if (!file.exists(path)) return(NULL)
    d <- read.table(path, header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","M","U","ctx","tri"))
    d <- d[d$ctx=="CG" & (d$M+d$U)>=4,]
    d$beta <- d$M/(d$M+d$U)
    GRanges(d$chr, IRanges(d$pos,d$pos), beta=d$beta)
  })
  grs <- Filter(Negate(is.null), grs)
  if (length(grs)==0) return(NULL)
  do.call(c, grs)
}

plots <- list()
for (ct in CONTRASTS) {
  message("Processing: ", ct$name)
  gr <- readRDS(ct$rds)
  gr <- gr[gr$context=="CG" & as.character(seqnames(gr))=="chr5"]
  if (length(gr)==0) { message("  no chr5 DMRs"); next }
  # subsample to avoid OOM on large contrasts
  if (length(gr)>5000) { set.seed(42); gr <- gr[sample(length(gr),5000)] }
  message("  using ", length(gr), " DMRs")
  cpg <- load_cx_chr5(ct$cond)
  if (is.null(cpg)) { message("  no CpG data"); next }
  ov <- findOverlaps(cpg, gr)
  dmr_beta <- tapply(cpg$beta[queryHits(ov)], subjectHits(ov), mean, na.rm=TRUE)
  dmr_df <- data.frame(
    beta=as.numeric(dmr_beta),
    direction=ifelse(gr$regionType[as.integer(names(dmr_beta))]=="gain","Hyper","Hypo"))
  set.seed(42)
  bg_idx <- sample(seq_along(cpg), min(20000, length(cpg)))
  bg_df <- data.frame(beta=cpg$beta[bg_idx], direction="Background")
  plot_df <- rbind(dmr_df, bg_df)
  plot_df$direction <- factor(plot_df$direction, levels=c("Background","Hypo","Hyper"))
  p <- ggplot(plot_df, aes(x=beta, fill=direction)) +
    geom_density(alpha=0.6, colour="grey25") +
    scale_fill_manual(values=c("Background"="#BDBDBD","Hypo"="#1F3A5F","Hyper"="#C0392B"), name=NULL) +
    scale_x_continuous(labels=scales::percent_format(1), limits=c(0,1)) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold")) +
    labs(title=paste0("Methylation distribution: ", ct$name, " (chr5)"),
         subtitle="Mean beta per DMR vs background CpGs",
         x="CpG methylation", y="Density")
  ggsave(file.path(OUT_DIR, paste0("meth_dist_",ct$name,".pdf")),
         p, width=7, height=4.5, device=cairo_pdf)
  plots[[ct$name]] <- p
  message("  saved: meth_dist_", ct$name)
}

if (length(plots)>=2) {
  library(patchwork)
  combined <- wrap_plots(plots, ncol=1) +
    plot_annotation(title="Methylation level distributions at DMR loci",
                    theme=theme(plot.title=element_text(face="bold",size=13)))
  ggsave(file.path(OUT_DIR,"meth_dist_combined.pdf"),
         combined, width=7, height=4.5*length(plots), device=cairo_pdf)
}
message("Done. Outputs in: ", OUT_DIR)
