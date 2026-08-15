.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
source("scripts/pipeline/00_sma_palette.R")

OUT_DIR <- "results/global_methylation"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ---- Safety check: does any condition show a massive, unbounded global
# methylation shift compared to the others, beyond the expected VPA effect?
# Computed per REPLICATE (not pooled) so we get real sample-level variance,
# which is what a safety/global-shift plot should show.

BY_CHR_DIR <- "results/alignments/bs/by_chr"
CONDITIONS <- c("ASO_CTRL", "Scramble_CTRL", "ASO_VPA", "Scramble_VPA")
REPS <- 1:3
CHROMS <- paste0("chr", 1:22)  # autosomes only, consistent with rest of pipeline

COND_COLOURS <- c(
  ASO_CTRL      = "#1B4F8A",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500"
)

compute_sample_mean <- function(condition, rep) {
  total_M <- 0
  total_N <- 0
  for (ch in CHROMS) {
    f <- file.path(BY_CHR_DIR, sprintf("%s_%d_%s.CpG_report.txt.gz", condition, rep, ch))
    if (!file.exists(f)) next
    d <- read.table(gzfile(f), header=FALSE, sep="\t",
                    col.names=c("chr","pos","strand","countM","countU","context","tri"),
                    colClasses=c("character","integer","character","integer",
                                 "integer","character","character"))
    d <- d[d$context=="CG", ]
    d <- d[d$countM + d$countU >= 4, ]
    total_M <- total_M + sum(d$countM)
    total_N <- total_N + sum(d$countM + d$countU)
  }
  if (total_N == 0) return(NA)
  total_M / total_N
}

message("Computing genome-wide mean methylation per sample...")
results <- data.frame()
for (cond in CONDITIONS) {
  for (r in REPS) {
    message("  ", cond, " rep ", r)
    m <- compute_sample_mean(cond, r)
    results <- rbind(results, data.frame(condition=cond, replicate=r, mean_methylation=m))
    message(sprintf("    mean methylation = %.4f", m))
  }
}

write.csv(results, file.path(OUT_DIR, "global_methylation_per_sample.csv"), row.names=FALSE)
message("\nSaved: global_methylation_per_sample.csv")

summary_stats <- aggregate(mean_methylation ~ condition, data=results, FUN=function(x) c(mean=mean(x), sd=sd(x)))
summary_stats <- do.call(data.frame, summary_stats)
names(summary_stats) <- c("condition", "mean", "sd")
summary_stats$condition <- factor(summary_stats$condition, levels=CONDITIONS)
results$condition <- factor(results$condition, levels=CONDITIONS)

p <- ggplot(summary_stats, aes(x=condition, y=mean, colour=condition)) +
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd), size=1, linewidth=1.1) +
  geom_jitter(data=results, aes(x=condition, y=mean_methylation), width=0.06, size=1.8, alpha=0.5, inherit.aes=FALSE, colour="grey30") +
  scale_colour_manual(values=COND_COLOURS, guide="none") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=14),
        plot.subtitle=element_text(size=10.5, colour="grey30"),
        axis.text.x=element_text(size=11)) +
  labs(title="Genome-wide mean CpG methylation by condition",
       subtitle="n=3 biological replicates per condition, autosomes, sites with >=4x coverage\nVPA-treated groups show lower genome-wide mean CpG methylation than controls; ASO_CTRL is close to Scramble_CTRL",
       x=NULL, y="Mean CpG methylation (genome-wide)")

ggsave(file.path(OUT_DIR, "global_methylation_safety_boxplot.pdf"), p,
       width=9, height=6.5, device=cairo_pdf)
message("saved: global_methylation_safety_boxplot.pdf")

cat("\n=== Summary ===\n")
print(aggregate(mean_methylation ~ condition, data=results,
                FUN=function(x) c(mean=mean(x), sd=sd(x))))
