#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
})

# Additional QC plots:
# 1. M-bias (positional methylation bias along reads)
# 2. Duplication rates
# 3. Bisulfite conversion efficiency (from CHH context)
# 4. Per-condition methylation violin at high-confidence DMR loci

BS_DIR <- "results/alignments/bs"
OUT    <- "results/qc/additional_qc"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

SAMPLES <- c(
  "ASO_CTRL_1","ASO_CTRL_2","ASO_CTRL_3",
  "ASO_VPA_1","ASO_VPA_2","ASO_VPA_3",
  "Scramble_CTRL_1","Scramble_CTRL_2","Scramble_CTRL_3",
  "Scramble_VPA_1","Scramble_VPA_2","Scramble_VPA_3"
)
CONDITION <- setNames(
  c(rep("ASO_CTRL",3), rep("ASO_VPA",3),
    rep("Scramble_CTRL",3), rep("Scramble_VPA",3)),
  SAMPLES)

GROUP_COLS <- c(
  ASO_CTRL="#1B4F8A", ASO_VPA="#B2182B",
  Scramble_VPA="#F0A500", Scramble_CTRL="#6B7280")

# M-bias
message("reading M-bias files...")
read_mbias <- function(sample) {
  f     <- file.path(BS_DIR, paste0(sample, "_bismark.deduplicated.M-bias.txt"))
  lines <- readLines(f)
  ctx_starts <- grep("^(CpG|CHG|CHH) context", lines)
  dfs <- list()
  for (i in seq_along(ctx_starts)) {
    end        <- if (i < length(ctx_starts)) ctx_starts[i+1] else length(lines)+1
    data_lines <- lines[(ctx_starts[i]+3):(end-1)]
    data_lines <- data_lines[nchar(trimws(data_lines)) > 0]
    df <- read.table(text=data_lines, header=FALSE,
                     col.names=c("position","methylated","unmethylated",
                                 "pct_methylation","coverage"))
    df$context <- gsub(" context \\(.*\\)", "", lines[ctx_starts[i]])
    df$read    <- gsub(".*\\((.*)\\)", "\\1", lines[ctx_starts[i]])
    df$sample  <- sample
    dfs[[i]]   <- df
  }
  do.call(rbind, dfs)
}

mbias_all <- do.call(rbind, lapply(SAMPLES, function(s) {
  message("  ", s); read_mbias(s) }))
mbias_all$condition <- CONDITION[mbias_all$sample]
mbias_cpg <- mbias_all[mbias_all$context == "CpG", ]

p_mbias <- ggplot(mbias_cpg, aes(x=position, y=pct_methylation,
                                  colour=condition, group=sample)) +
  geom_line(alpha=0.7, linewidth=0.6) +
  facet_wrap(~read, ncol=1) +
  scale_colour_manual(values=GROUP_COLS) +
  labs(title="M-bias — CpG context",
       x="Position in read (bp)", y="% CpG methylation") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"), strip.text=element_text(face="bold"))
ggsave(file.path(OUT, "mbias_CpG_all_samples.pdf"), p_mbias, width=10, height=8)
message("saved M-bias plot")

# duplication rates
message("reading duplication rates...")
dup_df <- do.call(rbind, lapply(SAMPLES, function(s) {
  f        <- file.path(BS_DIR, paste0(s, "_bismark.deduplication_report.txt"))
  pct_line <- grep("Total number duplicated", readLines(f), value=TRUE)
  pct      <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", pct_line))
  data.frame(sample=s, condition=CONDITION[s], dup_rate=pct)
}))
dup_df$sample <- factor(dup_df$sample, levels=SAMPLES)

p_dup <- ggplot(dup_df, aes(x=sample, y=dup_rate, fill=condition)) +
  geom_col() +
  geom_text(aes(label=paste0(round(dup_rate,1),"%")), vjust=-0.3, size=3) +
  scale_fill_manual(values=GROUP_COLS) +
  scale_y_continuous(limits=c(0, max(dup_df$dup_rate)*1.2)) +
  labs(title="PCR duplication rates", x=NULL, y="Duplication rate (%)") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        plot.title=element_text(face="bold"))
ggsave(file.path(OUT, "duplication_rates.pdf"), p_dup, width=10, height=5)
message("saved duplication rates")

# bisulfite conversion efficiency from CHH context
message("reading conversion efficiency...")
conv_df <- do.call(rbind, lapply(SAMPLES, function(s) {
  f     <- file.path(BS_DIR, paste0(s, "_bismark.deduplicated_splitting_report.txt"))
  lines <- readLines(f)
  chh_m <- as.numeric(gsub(".*:\t", "",
             grep("Total methylated C's in CHH", lines, value=TRUE)))
  chh_u <- as.numeric(gsub(".*:\t", "",
             grep("Total C to T conversions in CHH", lines, value=TRUE)))
  data.frame(sample=s, condition=CONDITION[s],
             conversion_efficiency_pct=round(chh_u*100/(chh_u+chh_m), 1))
}))
conv_df$sample <- factor(conv_df$sample, levels=SAMPLES)
write.csv(conv_df, file.path(OUT, "bisulfite_conversion_efficiency.csv"), row.names=FALSE)

p_conv <- ggplot(conv_df, aes(x=sample, y=conversion_efficiency_pct, fill=condition)) +
  geom_col() +
  geom_text(aes(label=paste0(conversion_efficiency_pct,"%")), vjust=-0.3, size=3) +
  scale_fill_manual(values=GROUP_COLS) +
  scale_y_continuous(limits=c(94,100)) +
  labs(title="Bisulfite conversion efficiency",
       subtitle="estimated from CHH context (expected ~0% true methylation in human somatic cells)",
       x=NULL, y="Conversion efficiency (%)") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        plot.title=element_text(face="bold"))
ggsave(file.path(OUT, "bisulfite_conversion_efficiency.pdf"), p_conv, width=10, height=5)
message("saved conversion efficiency")

# per-condition methylation violin at high-confidence DMR loci
message("generating methylation violin...")
get_meth <- function(contrast, col) {
  d <- readRDS(paste0("results/dmr/dmr_", contrast, ".rds"))
  d <- d[d$cytosinesCount >= 6]
  data.frame(methylation=mcols(d)[[col]], condition=col)
}

violin_df <- rbind(
  data.frame(methylation=local({
    d <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
    d[d$cytosinesCount>=6]$proportion1}), condition="ASO_CTRL"),
  data.frame(methylation=local({
    d <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
    d[d$cytosinesCount>=6]$proportion2}), condition="Scramble_CTRL"),
  data.frame(methylation=local({
    d <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
    d[d$cytosinesCount>=6]$proportion1}), condition="ASO_VPA"),
  data.frame(methylation=local({
    d <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
    d[d$cytosinesCount>=6]$proportion1}), condition="Scramble_VPA")
)
violin_df$condition <- factor(violin_df$condition,
  levels=c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA"))

p_violin <- ggplot(violin_df, aes(x=condition, y=methylation, fill=condition)) +
  geom_violin(trim=FALSE, alpha=0.8) +
  geom_boxplot(width=0.08, fill="white", outlier.size=0.2, outlier.alpha=0.2) +
  scale_fill_manual(values=GROUP_COLS) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(title="CpG methylation at high-confidence DMR loci",
       x=NULL, y="CpG methylation proportion") +
  theme_classic(base_size=12) +
  theme(legend.position="none", plot.title=element_text(face="bold"),
        axis.text.x=element_text(size=11, face="bold"))
ggsave(file.path(OUT, "per_condition_methylation_violin.pdf"), p_violin, width=8, height=6)
message("saved violin")
message("\ndone. outputs in: ", OUT)
