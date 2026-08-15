.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/smn2_enhancer"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

SMN2_START <- 70049638; SMN2_END <- 70078522
REGION_S   <- 69000000; REGION_E  <- 72000000
MIN_VIS    <- 5000
EFZ_START  <- 69530000; EFZ_END <- 71010000  # 1.48 Mb enhancer-free zone

CONTRASTS <- c(
  "ASO_CTRL_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_VPA_vs_Scramble_CTRL",
  "ASO_VPA_vs_ASO_CTRL",
  "ASO_VPA_vs_Scramble_VPA"
)

# ---- load enhancers ----
message("loading H9 enhancers...")
bed <- read.table(pipe("zcat data/reference/H9_predicted_non_promoter_non_fragments.bed.gz"),
  header=TRUE, sep="\t", stringsAsFactors=FALSE)
bed5 <- bed[bed$seqnames=="chr5" & bed$start>=REGION_S & bed$end<=REGION_E,]
enh_gr <- GRanges(bed5$seqnames, IRanges(bed5$start, bed5$end))
enh_df <- data.frame(start=start(enh_gr), end=end(enh_gr))
message("  enhancers in region: ", length(enh_gr))

smn2_df <- data.frame(start=SMN2_START, end=SMN2_END)

# ---- load DMRs ----
dmr_list <- list()
for (ct in CONTRASTS) {
  rds <- paste0("results/dmr/dmr_",ct,".rds")
  if (!file.exists(rds)) next
  gr <- readRDS(rds)
  gr <- gr[as.character(seqnames(gr))=="chr5" &
           start(gr)>=REGION_S & end(gr)<=REGION_E &
           gr$context=="CG"]
  if (length(gr)==0) {
    message(sprintf("  %-44s 0 DMRs in window", ct))
    next
  }
  if (length(gr)>50000) { set.seed(42); gr <- gr[sample(length(gr),50000)] }
  dmr_list[[ct]] <- data.frame(
    start=start(gr), end=end(gr),
    direction=ifelse(gr$regionType=="gain","hypo","hyper"),
    contrast=ct, stringsAsFactors=FALSE)
  message(sprintf("  %-44s %6d hypo  %6d hyper", ct,
    sum(gr$regionType=="gain"), sum(gr$regionType=="loss")))
}
dmr_df <- do.call(rbind, dmr_list)

cat("\nhyper/hypo counts per contrast in window:\n")
print(table(dmr_df$contrast, dmr_df$direction, useNA="ifany"))

# pad + y positions
dmr_plot <- dmr_df %>%
  mutate(
    contrast  = factor(contrast, levels=CONTRASTS),
    y_centre  = as.numeric(contrast),
    pad       = pmax(0, (MIN_VIS-(end-start))/2),
    xmin_v    = start - pad,
    xmax_v    = end   + pad,
    y_lo      = ifelse(direction=="hypo",  y_centre-0.40, y_centre),
    y_hi      = ifelse(direction=="hypo",  y_centre,       y_centre+0.40)
  )

x_lims  <- c(REGION_S, REGION_E)
x_scale <- scale_x_continuous(
  limits=x_lims,
  labels=function(x) sprintf("%.1f Mb", x/1e6),
  expand=expansion(mult=c(0.01,0.03)))

# ---- Panel A: enhancer track ----
p_enh <- ggplot() +
  geom_rect(data=enh_df,
    aes(xmin=pmax(start-1000,x_lims[1]), xmax=pmin(end+1000,x_lims[2]),
        ymin=0.3, ymax=0.7),
    fill="#2E7D32", colour=NA) +
  geom_rect(data=smn2_df,
    aes(xmin=start, xmax=end, ymin=0.1, ymax=0.9),
    fill="#C0392B", colour=NA) +
  annotate("text", x=(SMN2_START+SMN2_END)/2, y=1.15,
    label="SMN2", colour="#C0392B", fontface="bold", size=4.5) +
  annotate("segment", x=EFZ_START, xend=EFZ_END, y=-0.2, yend=-0.2,
    arrow=arrow(ends="both", length=unit(0.15,"cm")), colour="grey40") +
  annotate("text", x=(EFZ_START+EFZ_END)/2, y=-0.55,
    label="1.48 Mb enhancer-free zone", colour="grey25", size=3.5) +
  scale_y_continuous(limits=c(-0.9,1.4), expand=c(0,0)) +
  x_scale +
  labs(y=sprintf("H9 enhancers\n(n=%d)", length(enh_gr))) +
  theme_classic(base_size=12) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=9, colour="grey25", angle=0, vjust=0.5),
        plot.margin=margin(5,5,0,5))

# ---- Panel B: DMR tracks ----
p_dmr <- ggplot() +
  geom_blank(data=data.frame(
    contrast=factor(CONTRASTS, levels=CONTRASTS),
    y=seq_along(CONTRASTS)),
    aes(y=y)) +
  geom_hline(data=data.frame(y=(seq_along(CONTRASTS)-0.5)[-1]),
    aes(yintercept=y), colour="grey90", linewidth=0.3) +
  geom_rect(data=filter(dmr_plot, direction=="hypo"),
    aes(xmin=xmin_v, xmax=xmax_v, ymin=y_lo, ymax=y_hi, fill="Hypo"),
    colour=NA) +
  geom_rect(data=filter(dmr_plot, direction=="hyper"),
    aes(xmin=xmin_v, xmax=xmax_v, ymin=y_lo, ymax=y_hi, fill="Hyper"),
    colour=NA) +
  scale_fill_manual(name="DMR direction",
    values=c("Hyper"="#C0392B","Hypo"="#1F3A5F")) +
  scale_y_continuous(
    breaks=seq_along(CONTRASTS), labels=CONTRASTS,
    limits=c(0.5, length(CONTRASTS)+0.5), expand=c(0,0)) +
  x_scale +
  labs(x="chr5 position", y=NULL) +
  theme_classic(base_size=12) +
  theme(panel.grid.major.x=element_line(colour="grey92"),
        axis.text.y=element_text(face="bold", size=9),
        legend.position="bottom",
        plot.margin=margin(0,5,5,5))

# ---- combine ----
final <- p_enh / p_dmr +
  plot_layout(heights=c(0.6,4)) +
  plot_annotation(
    title="H9 enhancers and DMRs around SMN2 (chr5:69-72 Mb)",
    subtitle="Each row = one contrast; upper half=hyper (red), lower half=hypo (blue); ticks padded to 5 kb",
    theme=theme(plot.title=element_text(face="bold",size=13)))

ggsave(file.path(OUT_DIR,"SMN2_enhancer_map.pdf"),
       final, width=12, height=8, device=cairo_pdf)
message("saved: SMN2_enhancer_map.pdf")
message("done. outputs in: ", OUT_DIR)
