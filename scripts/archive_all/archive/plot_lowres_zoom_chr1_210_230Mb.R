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

# DMR-dense region: chr1:210-230 Mb (165 DMRs in 220-230Mb bin alone)
REGION_START <- 210e6
REGION_END   <- 230e6
BIN_SIZE     <- 50000

message("Reading chr1:210-230 Mb for all samples...")
all_bins <- mapply(function(s, g) {
  f <- file.path(cov_dir, paste0(s, "_chr1.CpG_report.txt.gz"))
  message("  Reading: ", s)
  df <- read.table(gzfile(f), header=FALSE, sep="\t",
                   colClasses=c("NULL","integer","NULL",
                                "integer","integer","NULL","NULL"))
  colnames(df) <- c("pos","M","U")
  # Filter to region only before any further processing
  df <- df[df$pos >= REGION_START & df$pos <= REGION_END, ]
  cov <- df$M + df$U
  df <- df[cov >= 1, ]
  if (nrow(df) == 0) return(NULL)
  df$bin <- floor(df$pos / BIN_SIZE) * BIN_SIZE
  agg <- aggregate(cbind(M, U) ~ bin, data=df, FUN=sum)
  agg$meth   <- agg$M / (agg$M + agg$U)
  agg$sample <- s
  agg$group  <- g
  agg
}, samples, groups, SIMPLIFY=FALSE)

dat <- do.call(rbind, all_bins)
group_means <- aggregate(meth ~ bin + group, data=dat, FUN=mean)

# Group colours — consistent with full chr1 plots
GROUP_COLS <- c(
  "ASO_CTRL"      = "#2E9B6F",
  "ASO_VPA"       = "#D94F3D",
  "Scramble_CTRL" = "#1D6FA4",
  "Scramble_VPA"  = "#F0A500"
)

contrasts <- list(
  list(groups  = c("ASO_VPA","ASO_CTRL"),
       tag     = "ASO_effect"),
  list(groups  = c("Scramble_VPA","Scramble_CTRL"),
       tag     = "VPA_effect"),
  list(groups  = c("ASO_VPA","Scramble_CTRL"),
       tag     = "combined_effect")
)

for (ct in contrasts) {
  message("Plotting: ", ct$tag)
  d_rep  <- dat[dat$group %in% ct$groups, ]
  d_mean <- group_means[group_means$group %in% ct$groups, ]
  d_mean$group <- factor(d_mean$group, levels = ct$groups)
  d_rep$group  <- factor(d_rep$group,  levels = ct$groups)
  cols <- GROUP_COLS[ct$groups]

  p <- ggplot() +
    geom_line(data = d_rep,
              aes(x = bin/1e6, y = meth, group = sample, colour = group),
              linewidth = 0.35, alpha = 0.4) +
    geom_line(data = d_mean,
              aes(x = bin/1e6, y = meth, colour = group, group = group),
              linewidth = 1.3) +
    scale_colour_manual(values = cols, name = "Group") +
    scale_x_continuous(breaks = seq(210, 230, 2),
                       labels = paste0(seq(210, 230, 2), " Mb")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title    = paste0("CpG methylation chr1:210-230 Mb (50 kb bins)\n",
                        gsub("_", " ", ct$tag), ": ",
                        paste(ct$groups, collapse = " vs ")),
      subtitle = paste0("DMR-dense region (165 DMRs in 220-230 Mb). ",
                        "Thick lines = group means; thin = replicates."),
      x = "chr1 position (Mb)",
      y = "Mean CpG methylation"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position  = "bottom",
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )

  outfile <- file.path(out_dir,
    paste0("lowres_zoom_chr1_210_230Mb_", ct$tag, ".pdf"))
  ggsave(outfile, p, width = 10, height = 4.5)
  message("  Saved: ", outfile)
}
message("Done.")
