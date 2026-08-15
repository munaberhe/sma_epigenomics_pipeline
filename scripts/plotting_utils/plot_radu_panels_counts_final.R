.libPaths("~/R/library")
library(ggplot2)
library(patchwork)
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

df <- read.csv("results/dmr_benchmark_mac/benchmark_summary_wide.csv")
df <- df[df$mode=="strict",]
df$method <- factor(df$method, levels=c("Bins","Neighbourhood","Noise_filter"),
                    labels=c("DMRcaller-B","DMRcaller-NB","DMRcaller-NF"))
df$scramble_method <- factor(df$scramble_method,
  levels=c("Read count permutation","Label-swap","Stratified scramble"))
df$window_size_num <- as.numeric(as.character(df$window_size_num))

# exclude NB from coverage panels
df_plot <- df[df$method != "DMRcaller-NB",]

NAVY      <- "#1F3A5F"; NAVY_L <- "#7FB0D3"
RED       <- "#C0392B"; RED_L  <- "#E8A39B"

make_panels <- function(scramble_label, tag) {
  d <- df_plot[df_plot$scramble_method == scramble_label,]
  if (nrow(d)==0) { message("skip: ", scramble_label); return(NULL) }

  # Panel A: DMR counts
  d_long <- rbind(
    data.frame(method=d$method, x=d$window_size_num*0.95, y=d$n_real,      series="real"),
    data.frame(method=d$method, x=d$window_size_num*1.05, y=d$n_scrambled, series="scrambled")
  )
  d_long$key <- paste0(d_long$method, " (", d_long$series, ")")
  key_cols <- c("DMRcaller-B (real)"=NAVY, "DMRcaller-B (scrambled)"=NAVY_L,
                "DMRcaller-NF (real)"=RED,  "DMRcaller-NF (scrambled)"=RED_L)
  key_lty  <- c("DMRcaller-B (real)"="solid","DMRcaller-B (scrambled)"="longdash",
                "DMRcaller-NF (real)"="solid","DMRcaller-NF (scrambled)"="longdash")

  pA <- ggplot(d_long, aes(x=x, y=y, colour=key, linetype=key, group=key)) +
    geom_line(linewidth=1) + geom_point(size=3) +
    scale_colour_manual(values=key_cols, name=NULL) +
    scale_linetype_manual(values=key_lty, name=NULL) +
    scale_x_log10(breaks=c(100,200,300,500,1000,2000)) +
    scale_y_continuous(labels=scales::comma) +
    theme_minimal(base_size=12) +
    theme(plot.title=element_text(face="bold"), legend.position=c(0.02,0.98),
          legend.justification=c(0,1), panel.grid.minor=element_blank()) +
    labs(title=paste0("A  Number of DMRs (strict) — ", scramble_label),
         x="Bin/window size (bp)", y="Number of DMRs")

  # Panel B: delta counts
  pB <- ggplot(d, aes(x=window_size_num, y=n_real-n_scrambled,
                      colour=method, group=method)) +
    geom_hline(yintercept=0, linetype="dotted", colour="grey60") +
    geom_line(linewidth=1) + geom_point(size=3) +
    scale_colour_manual(values=c("DMRcaller-B"=NAVY,"DMRcaller-NF"=RED), name=NULL) +
    scale_x_log10(breaks=c(100,200,300,500,1000,2000)) +
    scale_y_continuous(labels=scales::comma) +
    theme_minimal(base_size=12) +
    theme(plot.title=element_text(face="bold"), legend.position=c(0.98,0.98),
          legend.justification=c(1,1), panel.grid.minor=element_blank()) +
    labs(title="B  Delta DMRs (real minus scrambled, strict)",
         x="Bin/window size (bp)", y="Delta DMRs (real − scrambled)")

  combined <- (pA | pB) +
    plot_annotation(caption=paste0("null model: ", scramble_label,
      " | DMRcaller-NB excluded (count invariant to window size)"),
      theme=theme(plot.caption=element_text(colour="grey40", face="italic")))

  out <- file.path("results/dmr_benchmark/plots",
                   paste0("radu_panels_counts_", tag, ".pdf"))
  dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
  ggsave(out, combined, width=14, height=5.5, device=cairo_pdf)
  message("Saved: ", basename(out))
}

make_panels("Read count permutation", "read_count_permutation")
make_panels("Label-swap",             "label_swap")
make_panels("Stratified scramble",    "stratified_scramble")
message("Done.")
