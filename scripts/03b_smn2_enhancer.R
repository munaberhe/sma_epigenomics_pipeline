# SMN2 locus enhancer analysis
# H9 ESC predicted non-promoter enhancers (Khizra) vs ASO/VPA DMRs.
# question: do any of our DMRs overlap predicted enhancers at SMN2?
# answer: SMN2 sits in a ~1.48 Mb enhancer desert (chr5:69.53-71.01 Mb).
# a handful of methylation gains in VPA contrasts fall on enhancers ~1 Mb
# upstream of SMN2, not at SMN2 itself.
.libPaths("~/R/library")
suppressPackageStartupMessages({ library(GenomicRanges); library(ggplot2) })
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/smn2_enhancer"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
SMN2_START <- 70049638; SMN2_END <- 70078522
REGION_S <- 69000000;   REGION_E <- 72000000
CONTRASTS <- c("ASO_CTRL_vs_Scramble_CTRL","ASO_VPA_vs_Scramble_CTRL",
               "ASO_VPA_vs_ASO_CTRL","ASO_VPA_vs_Scramble_VPA",
               "Scramble_VPA_vs_Scramble_CTRL")
band_of <- function(s, e) {
  d <- pmin(abs(s - SMN2_START), abs(e - SMN2_END))
  inside <- (e >= SMN2_START) & (s <= SMN2_END); d[inside] <- 0
  cut(d, breaks=c(-1,0,100000,500000,1000000,Inf),
      labels=c("at SMN2","<100 kb","100-500 kb","500 kb-1 Mb",">1 Mb"))
}
message("loading H9 enhancers...")
bed <- read.table(pipe("zcat data/reference/H9_predicted_non_promoter_non_fragments.bed.gz"),
                  header=TRUE, sep="\t", stringsAsFactors=FALSE)
bed5 <- bed[bed$seqnames=="chr5" & bed$start>=REGION_S & bed$end<=REGION_E,]
enh_gr <- GRanges(bed5$seqnames, IRanges(bed5$start, bed5$end), type=bed5$type)
up_idx   <- which(end(enh_gr)   < SMN2_START)
down_idx <- which(start(enh_gr) > SMN2_END)
nearest_up   <- if(length(up_idx))   max(end(enh_gr)[up_idx])     else NA
nearest_down <- if(length(down_idx)) min(start(enh_gr)[down_idx]) else NA
message("  enhancers in region: ", length(enh_gr))
if(!is.na(nearest_up))
  message("  nearest upstream end: ", nearest_up,
          " (", round((SMN2_START-nearest_up)/1e3,1), " kb upstream)")
if(!is.na(nearest_down))
  message("  nearest downstream start: ", nearest_down,
          " (", round((nearest_down-SMN2_END)/1e3,1), " kb downstream)")
if(!is.na(nearest_up) && !is.na(nearest_down))
  message("  enhancer-free zone: ", round((nearest_down-nearest_up)/1e6,2), " Mb")
rows <- list(); summary_rows <- list()
for (ct in CONTRASTS) {
  dmr <- readRDS(paste0("results/dmr/dmr_",ct,".rds"))
  dmr5 <- dmr[as.character(seqnames(dmr))=="chr5"]
  dmr5 <- dmr5[start(dmr5)>=REGION_S & end(dmr5)<=REGION_E]
  key  <- paste(seqnames(dmr5), start(dmr5), end(dmr5), sep="_")
  dmr5 <- dmr5[!duplicated(key)]
  band_all <- band_of(start(dmr5), end(dmr5))
  hits <- subsetByOverlaps(dmr5, enh_gr)
  hkey <- paste(seqnames(hits), start(hits), end(hits), sep="_")
  hits <- hits[!duplicated(hkey)]
  message(sprintf("%-36s %3d DMRs, %d enhancer overlaps", ct, length(dmr5), length(hits)))
  if(length(hits)>0) {
    df <- as.data.frame(hits)
    df$contrast  <- ct
    df$methDiff  <- round(df$proportion2 - df$proportion1, 3)
    # gain = proportion2 > proportion1 = hypermethylation in cond2
    df$direction <- ifelse(df$methDiff > 0, "hyper", "hypo")
    df$distToSMN2_bp <- ifelse(df$end < SMN2_START, df$end - SMN2_START,
                          ifelse(df$start > SMN2_END, df$start - SMN2_END, 0))
    df$band <- as.character(band_of(df$start, df$end))
    rows[[ct]] <- df[,c("seqnames","start","end","regionType","pValue",
                        "proportion1","proportion2","methDiff","direction",
                        "distToSMN2_bp","band","contrast")]
  }
  hit_band <- if(length(hits)) band_of(start(hits), end(hits)) else factor(NULL)
  for(lvl in levels(band_all))
    summary_rows[[length(summary_rows)+1]] <- data.frame(
      contrast=ct, band=lvl,
      n_dmr=sum(band_all==lvl), n_enhancer_hit=sum(hit_band==lvl),
      stringsAsFactors=FALSE)
}
if(length(rows)>0) {
  hits_df <- do.call(rbind, rows)
  write.csv(hits_df, file.path(OUT_DIR,"SMN2_region_DMR_enhancer_overlaps.csv"), row.names=FALSE)
  message("saved: SMN2_region_DMR_enhancer_overlaps.csv (", nrow(hits_df), " rows)")
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_DIR,"SMN2_enhancer_overlap_by_band.csv"), row.names=FALSE)
print(summary_df, row.names=FALSE)
enh_df <- as.data.frame(enh_gr)
dmr_long <- list()
for(ct in CONTRASTS) {
  dmr <- readRDS(paste0("results/dmr/dmr_",ct,".rds"))
  dmr5 <- dmr[as.character(seqnames(dmr))=="chr5"]
  dmr5 <- dmr5[start(dmr5)>=REGION_S & end(dmr5)<=REGION_E]
  key <- paste(seqnames(dmr5), start(dmr5), end(dmr5), sep="_")
  dmr5 <- dmr5[!duplicated(key)]
  if(length(dmr5)==0) next
  dmr_long[[ct]] <- data.frame(start=start(dmr5), end=end(dmr5),
    methDiff=mcols(dmr5)$proportion2-mcols(dmr5)$proportion1, contrast=ct)
}
dmr_long <- if(length(dmr_long)) do.call(rbind, dmr_long) else NULL
p <- ggplot() +
  geom_rect(data=enh_df, aes(xmin=start/1e6, xmax=end/1e6, ymin=0.7, ymax=1.0),
            fill="#2E7D32", alpha=1.0, colour="#2E7D32", linewidth=0.3) +
  annotate("rect", xmin=SMN2_START/1e6, xmax=SMN2_END/1e6,
           ymin=0.6, ymax=1.1, fill="#C0392B", alpha=0.9) +
  annotate("text", x=(SMN2_START+SMN2_END)/2/1e6, y=1.18,
           label="SMN2", colour="#C0392B", size=3.5, fontface="bold")
if(!is.na(nearest_up) && !is.na(nearest_down)) {
  p <- p +
    annotate("segment", x=nearest_up/1e6, xend=nearest_down/1e6,
             y=0.45, yend=0.45, colour="grey40",
             arrow=arrow(ends="both", length=unit(0.12,"cm"))) +
    annotate("text", x=(nearest_up+nearest_down)/2/1e6, y=0.53,
             label=sprintf("enhancer-free zone: %.2f Mb",
                           (nearest_down-nearest_up)/1e6),
             size=3, colour="grey30")
}
if(!is.null(dmr_long) && nrow(dmr_long)>0) {
  dmr_long$direction <- ifelse(dmr_long$methDiff>0, "hyper", "hypo")
  dmr_long$y <- 0.3 - 0.05*(as.integer(factor(dmr_long$contrast, levels=CONTRASTS)))
  p <- p +
    geom_segment(data=dmr_long,
                 aes(x=start/1e6, xend=end/1e6, y=y, yend=y, colour=direction),
                 linewidth=2) +
    scale_colour_manual(values=c(hyper="#C0392B", hypo="#1F3A5F"), name="DMR")
  p <- p + geom_text(data=data.frame(contrast=CONTRASTS,
                       y=0.3-0.05*seq_along(CONTRASTS)),
                     aes(x=REGION_E/1e6+0.05, y=y, label=contrast),
                     hjust=0, size=2.5, colour="grey25")
}
p <- p +
  scale_x_continuous(limits=c(REGION_S/1e6, REGION_E/1e6+1.5),
                     breaks=seq(69,72,0.5),
                     labels=function(x) sprintf("%.1f Mb",x)) +
  scale_y_continuous(limits=c(-0.1,1.3)) +
  theme_classic(base_size=11) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        plot.title=element_text(face="bold"),
        plot.subtitle=element_text(colour="grey30",size=9),
        legend.position="bottom") +
  labs(title="H9 predicted enhancers and DMRs around SMN2 (chr5:69-72 Mb)",
       subtitle=sprintf("green=H9 enhancers (n=%d) | red=SMN2 | DMR tracks below (hyper=red, hypo=blue) | %.2f Mb enhancer-free zone",
         length(enh_gr),
         if(!is.na(nearest_up)&&!is.na(nearest_down)) (nearest_down-nearest_up)/1e6 else NA),
       x="chr5 position (Mb)", y=NULL)
ggsave(file.path(OUT_DIR,"SMN2_enhancer_map.pdf"), p,
       width=13, height=4.5, device=cairo_pdf)
message("saved: SMN2_enhancer_map.pdf")
message("done. outputs in: ", OUT_DIR)
