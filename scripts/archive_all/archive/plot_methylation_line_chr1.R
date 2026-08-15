#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages(library(ggplot2))

cov_dir <- "results/alignments/bs/by_chr"
out_dir <- "results/qc/for_meeting"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- c("ASO_CTRL_1","ASO_CTRL_2","ASO_CTRL_3",
             "ASO_VPA_1","ASO_VPA_2","ASO_VPA_3",
             "Scramble_CTRL_1","Scramble_CTRL_2","Scramble_CTRL_3",
             "Scramble_VPA_1","Scramble_VPA_2","Scramble_VPA_3")
groups <- c(rep("ASO_CTRL",3), rep("ASO_VPA",3),
            rep("Scramble_CTRL",3), rep("Scramble_VPA",3))

bin_size <- 1e6

message("Reading all samples...")
all_bins <- mapply(function(s, g) {
  f <- file.path(cov_dir, paste0(s, "_chr1.CpG_report.txt.gz"))
  message("  Reading: ", s)
  df <- read.table(gzfile(f), header=FALSE, sep="\t",
                   colClasses=c("NULL","integer","NULL",
                                "integer","integer","NULL","NULL"))
  colnames(df) <- c("pos","M","U")
  cov <- df$M + df$U
  df <- df[cov >= 1, ]
  df$bin <- floor(df$pos / bin_size) * bin_size
  agg <- aggregate(cbind(M, U) ~ bin, data=df, FUN=sum)
  agg$meth <- agg$M / (agg$M + agg$U)
  agg$sample <- s
  agg$group  <- g
  agg
}, samples, groups, SIMPLIFY=FALSE)

dat <- do.call(rbind, all_bins)
group_means <- aggregate(meth ~ bin + group, data=dat, FUN=mean)

contrasts <- list(
  list(groups=c("ASO_VPA","ASO_CTRL"),
       colours=c("ASO_VPA"="#D55E00","ASO_CTRL"="#0072B2"),
       title="Methylation profile along chr1 (1 Mb bins)\nASO_VPA vs ASO_CTRL",
       tag="ASO_effect"),
  list(groups=c("Scramble_VPA","Scramble_CTRL"),
       colours=c("Scramble_VPA"="#CC79A7","Scramble_CTRL"="#009E73"),
       title="Methylation profile along chr1 (1 Mb bins)\nScramble_VPA vs Scramble_CTRL",
       tag="VPA_effect"),
  list(groups=c("ASO_VPA","Scramble_CTRL"),
       colours=c("ASO_VPA"="#E69F00","Scramble_CTRL"="#56B4E9"),
       title="Methylation profile along chr1 (1 Mb bins)\nASO_VPA vs Scramble_CTRL",
       tag="combined_effect")
)

for (ct in contrasts) {
  message("Plotting: ", ct$tag)
  d_rep  <- dat[dat$group %in% ct$groups, ]
  d_mean <- group_means[group_means$group %in% ct$groups, ]
  d_mean$group <- factor(d_mean$group, levels=ct$groups)
  d_rep$group  <- factor(d_rep$group,  levels=ct$groups)

  p <- ggplot() +
    geom_line(data=d_rep,  aes(x=bin/1e6, y=meth, group=sample, colour=group),
              linewidth=0.3, alpha=0.4) +
    geom_line(data=d_mean, aes(x=bin/1e6, y=meth, colour=group, group=group),
              linewidth=1.2) +
    scale_colour_manual(values=ct$colours, name="Group") +
    scale_x_continuous(breaks=seq(0,250,50)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
    labs(title=ct$title,
         subtitle="Thick lines = group means; thin lines = individual replicates",
         x="chr1 position (Mb)", y="Mean CpG methylation") +
    theme_classic(base_size=12) +
    theme(legend.position="bottom",
          plot.title=element_text(face="bold", size=12),
          plot.subtitle=element_text(size=9, colour="grey40"))

  outfile <- file.path(out_dir, paste0("methylation_line_chr1_", ct$tag, ".pdf"))
  ggsave(outfile, p, width=10, height=4.5)
  message("  Saved: ", outfile)
}
message("Done.")
