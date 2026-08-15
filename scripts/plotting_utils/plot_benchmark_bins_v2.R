.libPaths("~/R/library")
library(ggplot2); library(patchwork)
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

df <- read.csv("results/dmr_benchmark_mac/parameter_benchmark_bins.csv")
df$window_size <- as.numeric(as.character(df$window_size))

NAVY <- "#1F3A5F"; RED <- "#C0392B"; GOLD <- "#D4A017"
mode_cols <- c(strict=NAVY, loose=adjustcolor(GOLD, 0.6))

# Panel A: strict mode only - real vs scrambled (loose is noise, don't show it)
d_strict <- df[df$mode=="strict",]
d_long <- rbind(
  data.frame(x=d_strict$window_size*0.95, y=d_strict$n_real,
             series="Real data", stringsAsFactors=FALSE),
  data.frame(x=d_strict$window_size*1.05, y=d_strict$n_scrambled,
             series="Label-swap null", stringsAsFactors=FALSE)
)

pA <- ggplot(d_long, aes(x=x, y=y, colour=series, linetype=series, group=series)) +
  geom_line(linewidth=1.2) + geom_point(size=3) +
  scale_colour_manual(values=c("Real data"=NAVY, "Label-swap null"=adjustcolor(NAVY,0.45)),
                      name=NULL) +
  scale_linetype_manual(values=c("Real data"="solid","Label-swap null"="dashed"),
                        name=NULL) +
  scale_x_log10(breaks=c(100,200,300,500,1000,2000)) +
  scale_y_log10(labels=scales::comma) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold", size=11),
        panel.grid.minor=element_blank(),
        legend.position="bottom") +
  labs(title="A  DMRcaller-B: DMR counts (strict mode, label-swap null)",
       x="Bin size (bp)", y="Number of DMRs (log scale)")

# Panel B: ratio for both modes
pB <- ggplot(df, aes(x=window_size, y=ratio, colour=mode, group=mode)) +
  geom_hline(yintercept=1, linetype="dotted", colour="grey50", linewidth=1) +
  geom_vline(xintercept=300, linetype="dashed", colour=RED, linewidth=0.9) +
  geom_line(linewidth=1.2) + geom_point(size=3) +
  annotate("text", x=380, y=max(df$ratio, na.rm=TRUE)*0.85,
           label="300 bp selected", colour=RED,
           hjust=0, size=3.2, fontface="bold") +
  annotate("text", x=380, y=1.05,
           label="null threshold", colour="grey50",
           hjust=0, size=2.8) +
  scale_colour_manual(values=c(strict=NAVY, loose=GOLD),
                      labels=c(strict="Strict mode",loose="Loose mode"),
                      name=NULL) +
  scale_x_log10(breaks=c(100,200,300,500,1000,2000)) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold", size=11),
        panel.grid.minor=element_blank(),
        legend.position="bottom") +
  labs(title="B  Signal-to-null ratio — strict crosses 1.0 at 300bp, loose remains below",
       x="Bin size (bp)", y="Real / Scrambled ratio")

combined <- (pA | pB) +
  plot_annotation(
    title="DMRcaller-B parameter selection: bins method, label-swap null",
    caption=paste0("Strict mode: real \u2265 null from 300 bp onward (ratio \u2265 1.0). ",
                   "Loose mode: real < null throughout (noise-dominated). 300 bp selected."),
    theme=theme(
      plot.title=element_text(face="bold", size=13),
      plot.caption=element_text(colour="grey40", face="italic", size=9,
                                hjust=0, margin=margin(t=8))))

ggsave("results/dmr_benchmark/plots/benchmark_bins_labelswap_v2.pdf",
       combined, width=13, height=6, device=cairo_pdf)
message("saved: benchmark_bins_labelswap_v2.pdf")
