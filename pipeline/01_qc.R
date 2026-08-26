#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

setwd("/data/home/bt25018/sma_epigenomics_pipeline")

BY_CHR  <- "results/alignments/bs/by_chr"
BS_DIR  <- "results/alignments/bs"
OUT_PCA <- "results/dmr_qc"
OUT_QC  <- "results/qc/additional_qc"
OUT_COV <- "results/qc/coverage_4lines"
MIN_COV <- 5
XMAX    <- 100

dir.create(OUT_PCA, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_QC,  recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_COV, recursive=TRUE, showWarnings=FALSE)

SAMPLES <- c(
  "ASO_CTRL_1","ASO_CTRL_2","ASO_CTRL_3",
  "ASO_VPA_1","ASO_VPA_2","ASO_VPA_3",
  "Scramble_CTRL_1","Scramble_CTRL_2","Scramble_CTRL_3",
  "Scramble_VPA_1","Scramble_VPA_2","Scramble_VPA_3"
)
CONDITIONS <- c("ASO_CTRL","ASO_VPA","Scramble_CTRL","Scramble_VPA")
CONDITION  <- c(rep("ASO_CTRL",3), rep("ASO_VPA",3),
                rep("Scramble_CTRL",3), rep("Scramble_VPA",3))
COND_NAMED <- setNames(CONDITION, SAMPLES)
REPS       <- 1:3
CHROMS     <- paste0("chr", c(1:22,"X","Y"))

GROUP_COLS <- c(
  ASO_CTRL="#1F3A5F", ASO_VPA="#C0392B",
  Scramble_VPA="#F0A500", Scramble_CTRL="#6B7280"
)
TRACK_COLS <- c(rep1="#7CB6D6", rep2="#4A8FB8", rep3="#1B5478", pooled="#C0392B")

message("reading chr1 for all 12 samples...")
meth_list <- list()
for (s in SAMPLES) {
  f <- file.path(BY_CHR, paste0(s, "_chr1.CpG_report.txt.gz"))
  gr <- readBismark(f)
  meth_list[[s]] <- gr[gr$context == "CG"]
  message("  ", s)
}

message("finding common CpGs (>= ", MIN_COV, "x)...")
common_pos <- Reduce(intersect, lapply(meth_list, function(gr) {
  gr <- gr[gr$readsN >= MIN_COV]
  paste0(as.character(seqnames(gr)), ":", start(gr))
}))
message("  common CpGs: ", length(common_pos))

message("building matrix...")
meth_mat <- do.call(cbind, lapply(SAMPLES, function(s) {
  gr  <- meth_list[[s]]
  pos <- paste0(as.character(seqnames(gr)), ":", start(gr))
  idx <- match(common_pos, pos)
  gr$readsM[idx] / gr$readsN[idx]
}))
colnames(meth_mat) <- SAMPLES
meth_mat <- meth_mat[complete.cases(meth_mat), ]
if (nrow(meth_mat) > 100000) {
  set.seed(42)
  meth_mat <- meth_mat[sample(nrow(meth_mat), 100000), ]
}

message("running PCA...")
pca    <- prcomp(t(meth_mat), scale.=FALSE, center=TRUE)
pca_df <- data.frame(
  PC1=pca$x[,1], PC2=pca$x[,2],
  sample=SAMPLES, condition=CONDITION,
  replicate=rep(c("1","2","3"), 4)
)
var_exp <- round(summary(pca)$importance[2,1:2]*100, 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour=condition,
                              shape=replicate, label=sample)) +
  geom_hline(yintercept=0, linetype="dotted", colour="grey60", linewidth=0.4) +
  geom_vline(xintercept=0, linetype="dotted", colour="grey60", linewidth=0.4) +
  geom_point(size=5, alpha=0.9,
             position=position_jitter(width=0.3, height=0.3, seed=42)) +
  geom_text_repel(size=3.2, colour="grey20", max.overlaps=20, box.padding=0.4,
                  segment.colour="grey70", segment.size=0.3,
                  position=position_jitter(width=0.3, height=0.3, seed=42)) +
  stat_ellipse(aes(group=condition), type="norm",
               linetype="dashed", linewidth=0.4, alpha=0.6) +
  scale_colour_manual(values=GROUP_COLS, name="Condition") +
  scale_shape_manual(values=c("1"=16,"2"=17,"3"=15), name="Replicate") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold"), legend.position="right") +
  labs(title="PCA of per-replicate CpG methylation (chr1)",
       x=paste0("PC1 (", var_exp[1], "%)"),
       y=paste0("PC2 (", var_exp[2], "%)"))

ggsave(file.path(OUT_PCA, "sample_PCA_12samples_chr1.pdf"),
       p_pca, width=9, height=7, device=cairo_pdf)
message("saved PCA")

violin_df <- data.frame(
  methylation=as.vector(meth_mat),
  sample=rep(colnames(meth_mat), each=nrow(meth_mat))
)
violin_df$condition <- CONDITION[match(violin_df$sample, SAMPLES)]
violin_df$sample    <- factor(violin_df$sample, levels=SAMPLES)
violin_df$condition <- factor(violin_df$condition,
  levels=c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA"))

p_violin <- ggplot(violin_df, aes(x=sample, y=methylation, fill=condition)) +
  geom_violin(trim=FALSE, alpha=0.8, linewidth=0.3) +
  geom_boxplot(width=0.07, fill="white",
               outlier.size=0.2, outlier.alpha=0.2) +
  scale_fill_manual(values=GROUP_COLS) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(title="Per-sample CpG methylation distribution (chr1)",
       x=NULL, y="CpG methylation proportion") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        plot.title=element_text(face="bold"), legend.position="top")

ggsave(file.path(OUT_PCA, "per_sample_methylation_violin_chr1.pdf"),
       p_violin, width=14, height=6, device=cairo_pdf)
message("saved violin")

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
mbias_all$condition <- COND_NAMED[mbias_all$sample]
mbias_cpg <- mbias_all[mbias_all$context == "CpG", ]

p_mbias <- ggplot(mbias_cpg, aes(x=position, y=pct_methylation,
                                   colour=condition, group=sample)) +
  geom_line(alpha=0.7, linewidth=0.6) +
  facet_wrap(~read, ncol=1) +
  scale_colour_manual(values=GROUP_COLS) +
  labs(title="M-bias - CpG context",
       x="Position in read (bp)", y="% CpG methylation") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"),
        strip.text=element_text(face="bold"))

ggsave(file.path(OUT_QC, "mbias_CpG_all_samples.pdf"),
       p_mbias, width=10, height=8)
message("saved M-bias")

message("reading duplication rates...")
dup_df <- do.call(rbind, lapply(SAMPLES, function(s) {
  f        <- file.path(BS_DIR, paste0(s, "_bismark.deduplication_report.txt"))
  pct_line <- grep("Total number duplicated", readLines(f), value=TRUE)
  pct      <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", pct_line))
  data.frame(sample=s, condition=COND_NAMED[s], dup_rate=pct)
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

ggsave(file.path(OUT_QC, "duplication_rates.pdf"), p_dup, width=10, height=5)
message("saved duplication rates")

message("reading conversion efficiency...")
conv_df <- do.call(rbind, lapply(SAMPLES, function(s) {
  f     <- file.path(BS_DIR,
             paste0(s, "_bismark.deduplicated_splitting_report.txt"))
  lines <- readLines(f)
  chh_m <- as.numeric(gsub(".*:\t", "",
             grep("Total methylated C's in CHH", lines, value=TRUE)))
  chh_u <- as.numeric(gsub(".*:\t", "",
             grep("Total C to T conversions in CHH", lines, value=TRUE)))
  data.frame(sample=s, condition=COND_NAMED[s],
             conversion_efficiency_pct=round(chh_u*100/(chh_u+chh_m), 1))
}))
conv_df$sample <- factor(conv_df$sample, levels=SAMPLES)
write.csv(conv_df,
          file.path(OUT_QC, "bisulfite_conversion_efficiency.csv"),
          row.names=FALSE)

p_conv <- ggplot(conv_df, aes(x=sample, y=conversion_efficiency_pct,
                               fill=condition)) +
  geom_col() +
  geom_text(aes(label=paste0(conversion_efficiency_pct,"%")),
            vjust=-0.3, size=3) +
  scale_fill_manual(values=GROUP_COLS) +
  scale_y_continuous(limits=c(94,100)) +
  labs(title="Bisulfite conversion efficiency", x=NULL,
       y="Conversion efficiency (%)") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        plot.title=element_text(face="bold"))

ggsave(file.path(OUT_QC, "bisulfite_conversion_efficiency.pdf"),
       p_conv, width=10, height=5)
message("saved conversion efficiency")

message("building coverage retention curves...")
read_cpg_report <- function(f) {
  d <- read.table(gzfile(f), header=FALSE, sep="\t",
    col.names=c("chr","pos","strand","countM","countU","context","tri"),
    colClasses=c("character","integer","character",
                 "integer","integer","character","character"))
  GRanges(seqnames=d$chr, ranges=IRanges(d$pos,d$pos), strand=d$strand,
          readsM=d$countM, readsN=d$countM+d$countU,
          context=d$context, trinucleotide_context=d$tri)
}

read_rep <- function(condition, rep, chroms) {
  files <- file.path(BY_CHR,
    sprintf("%s_%d_%s.CpG_report.txt.gz", condition, rep, chroms))
  files <- files[file.exists(files)]
  grs   <- lapply(files, read_cpg_report)
  do.call(c, unname(grs))
}

build_retention <- function(condition) {
  message("  ", condition)
  rep_data <- lapply(REPS, function(r) read_rep(condition, r, CHROMS))
  names(rep_data) <- paste0("rep", REPS)
  pooled   <- poolMethylationDatasets(GRangesList(rep_data))
  thresholds <- 0:XMAX
  do.call(rbind, lapply(c(paste0("rep",REPS),"pooled"), function(t) {
    cov <- if (t=="pooled") pooled$readsN else rep_data[[t]]$readsN
    data.frame(condition=condition, track=t, threshold=thresholds,
               retention=sapply(thresholds, function(k) mean(cov >= k)))
  }))
}

retention_all <- do.call(rbind, lapply(CONDITIONS, build_retention))
retention_all$condition <- factor(retention_all$condition, levels=CONDITIONS)
retention_all$track     <- factor(retention_all$track,
                                   levels=c("rep1","rep2","rep3","pooled"))

make_cov_panel <- function(cond) {
  ggplot(retention_all[retention_all$condition==cond,],
         aes(x=threshold, y=retention, colour=track)) +
    geom_line(linewidth=0.9) +
    coord_cartesian(xlim=c(0,XMAX), ylim=c(0,1)) +
    scale_colour_manual(values=TRACK_COLS) +
    geom_vline(xintercept=10, linetype=2, colour="grey40") +
    annotate("text", x=11, y=0.97, label="10x", hjust=0,
             colour="grey40", size=3.2) +
    labs(title=cond, x="Coverage threshold", y="Fraction CpGs >= threshold",
         colour=NULL) +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold", size=11),
          legend.position="right",
          panel.grid.major.y=element_line(colour="grey92"))
}

panels <- lapply(CONDITIONS, make_cov_panel)
combined_cov <- (panels[[1]] / panels[[2]] / panels[[3]] / panels[[4]]) +
  plot_layout(guides="collect") +
  plot_annotation(
    title="CpG coverage: replicates vs pooled",
    theme=theme(plot.title=element_text(face="bold", size=12))
  ) & theme(legend.position="right")

ggsave(file.path(OUT_COV, "coverage_4lines_per_condition.pdf"),
       combined_cov, width=8, height=12, device=cairo_pdf)
message("saved coverage plot")
message("done.")
