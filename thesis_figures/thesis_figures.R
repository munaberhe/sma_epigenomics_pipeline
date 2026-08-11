#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(reshape2)
  library(scales)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/thesis_figures"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

COND_COLS <- c(
  ASO_CTRL      = "#1F3A5F",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#C0392B",
  Scramble_VPA  = "#F0A500"
)
HYPO_COL  <- "#4393C3"
HYPER_COL <- "#D6604D"
CHROMS    <- paste0("chr", c(1:22, "X", "Y"))

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       label="ASO alone", short="ASO alone"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       label="VPA alone", short="VPA alone"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       label="ASO in VPA context", short="ASO in VPA"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       label="VPA in ASO context", short="VPA in ASO")
)

message("Loading DMR objects...")
dmrs <- lapply(CONTRASTS, function(ct) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) { message("  missing: ", rds); return(NULL) }
  readRDS(rds)
})
names(dmrs) <- sapply(CONTRASTS, `[[`, "name")

for (nm in names(dmrs)) {
  if (!is.null(dmrs[[nm]]))
    message("  ", nm, ": ", length(dmrs[[nm]]), " DMRs")
}

# Fig 5.1a — MDS plot (already exists, just copy)
file.copy("results/figures/qc_plots/mds_methylation.pdf",
          file.path(OUT, "Fig5.1a_MDS.pdf"), overwrite=TRUE)
message("Copied MDS plot")

# Fig 5.1b — Global CpG methylation violin per condition
message("Fig 5.1b: global methylation violin...")
meth_cache <- readRDS("results/dmr/meth_pooled_cache.rds")
violin_df <- do.call(rbind, lapply(names(meth_cache), function(cond) {
  m <- meth_cache[[cond]]
  m <- m[!is.na(m$readsN) & m$readsN >= 5]
  if (length(m) == 0) return(NULL)
  set.seed(42)
  idx <- sample(length(m), min(50000, length(m)))
  data.frame(
    methylation = m$readsM[idx] / m$readsN[idx],
    condition   = cond
  )
}))
violin_df$condition <- factor(violin_df$condition,
  levels=c("Scramble_CTRL","ASO_CTRL","Scramble_VPA","ASO_VPA"))

p_violin <- ggplot(violin_df, aes(x=condition, y=methylation, fill=condition)) +
  geom_violin(trim=FALSE, alpha=0.85, linewidth=0.3) +
  geom_boxplot(width=0.06, fill="white", outlier.size=0.3,
               outlier.alpha=0.2, linewidth=0.4) +
  scale_fill_manual(values=COND_COLS) +
  scale_y_continuous(labels=percent_format(accuracy=1),
                     breaks=seq(0,1,0.25)) +
  scale_x_discrete(labels=c(
    Scramble_CTRL="Scramble\nCTRL",
    ASO_CTRL="ASO\nCTRL",
    Scramble_VPA="Scramble\nVPA",
    ASO_VPA="ASO\nVPA")) +
  labs(title="Global CpG methylation distribution",
       subtitle="50,000 CpGs sampled per condition (≥5x coverage)",
       x=NULL, y="CpG methylation proportion") +
  theme_classic(base_size=11) +
  theme(legend.position="none",
        plot.title=element_text(face="bold"),
        axis.text.x=element_text(size=9))

ggsave(file.path(OUT, "Fig5.1b_global_methylation_violin.pdf"),
       p_violin, width=6, height=5, device="pdf")
message("Saved Fig5.1b")

# Fig 5.2 — DMR counts bar chart (Table 5.3 visual)
message("Fig 5.2: DMR counts...")
count_df <- do.call(rbind, lapply(seq_along(CONTRASTS), function(i) {
  ct <- CONTRASTS[[i]]
  d  <- dmrs[[ct$name]]
  if (is.null(d)) return(NULL)
  data.frame(
    contrast  = ct$short,
    direction = c("Hypomethylated","Hypermethylated"),
    count     = c(sum(mcols(d)$regionType=="gain"),
                  sum(mcols(d)$regionType=="loss")),
    signed    = c(sum(mcols(d)$regionType=="gain"),
                  -sum(mcols(d)$regionType=="loss"))
  )
}))
count_df$contrast <- factor(count_df$contrast,
  levels=c("ASO alone","VPA alone","ASO in VPA","VPA in ASO"))

p_counts <- ggplot(count_df, aes(x=contrast, y=signed, fill=direction)) +
  geom_bar(stat="identity", position="identity", width=0.6) +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey30") +
  scale_fill_manual(values=c(Hypomethylated=HYPO_COL,
                             Hypermethylated=HYPER_COL), name=NULL) +
  scale_y_continuous(labels=function(x) comma(abs(x))) +
  labs(title="DMR counts across four pairwise contrasts",
       x=NULL, y="Number of DMRs") +
  theme_classic(base_size=11) +
  theme(legend.position="top",
        plot.title=element_text(face="bold"),
        axis.text.x=element_text(size=9))

ggsave(file.path(OUT, "Fig5.2_DMR_counts.pdf"),
       p_counts, width=7, height=5, device="pdf")
message("Saved Fig5.2")

# Fig 5.3 — UpSet plot (already exists)
file.copy("results/figures/upset/upset_4contrasts.pdf",
          file.path(OUT, "Fig5.3_UpSet.pdf"), overwrite=TRUE)
message("Copied UpSet plot")

# Fig 5.4 — Circos (already exists)
file.copy("results/figures/genomic_distribution/circos_4contrasts_combined.pdf",
          file.path(OUT, "Fig5.4_Circos.pdf"), overwrite=TRUE)
message("Copied Circos plot")

# Fig 5.5 — Diverging chromosome bar (already exists)
file.copy("results/figures/genomic_distribution/chr_dmr_diverging_4contrasts.pdf",
          file.path(OUT, "Fig5.5_Chr_diverging.pdf"), overwrite=TRUE)
message("Copied diverging bar")

# Fig 5.6 — Methylation difference histograms (4 contrasts)
message("Fig 5.6: methylation diff histograms...")
hist_panels <- lapply(seq_along(CONTRASTS), function(i) {
  ct <- CONTRASTS[[i]]
  d  <- dmrs[[ct$name]]
  if (is.null(d)) return(NULL)
  df <- data.frame(
    diff      = mcols(d)$proportion1 - mcols(d)$proportion2,
    direction = ifelse(mcols(d)$regionType=="gain",
                       "Hypomethylated","Hypermethylated")
  )
  ggplot(df, aes(x=diff, fill=direction)) +
    geom_histogram(bins=50, colour="white", linewidth=0.1) +
    geom_vline(xintercept=0, linewidth=0.5, linetype="dashed") +
    scale_fill_manual(values=c(Hypomethylated=HYPO_COL,
                               Hypermethylated=HYPER_COL), name=NULL) +
    scale_x_continuous(limits=c(-1,1), labels=percent_format(1)) +
    scale_y_continuous(labels=comma) +
    labs(title=ct$label, x="methylation difference", y="DMRs") +
    theme_classic(base_size=9) +
    theme(plot.title=element_text(face="bold", size=9),
          legend.position="top")
})
fig56 <- (hist_panels[[1]] | hist_panels[[2]]) /
         (hist_panels[[3]] | hist_panels[[4]]) +
  plot_annotation(
    title="Methylation difference distributions across four pairwise contrasts",
    theme=theme(plot.title=element_text(face="bold", size=11))
  )
ggsave(file.path(OUT, "Fig5.6_meth_diff_histograms.pdf"),
       fig56, width=12, height=8, device="pdf")
message("Saved Fig5.6")

# Fig 5.7 — GO/KEGG (already exists)
file.copy("results/figures/gokegg_pairwise/GO_4contrasts_combined.pdf",
          file.path(OUT, "Fig5.7_GO.pdf"), overwrite=TRUE)
file.copy("results/figures/gokegg_pairwise/KEGG_4contrasts_combined.pdf",
          file.path(OUT, "Fig5.8_KEGG.pdf"), overwrite=TRUE)
message("Copied GO/KEGG plots")

# Fig 5.9 — Volcano (already exists)
file.copy("results/figures/volcano_plots/volcano_4contrasts.pdf",
          file.path(OUT, "Fig5.9_Volcano.pdf"), overwrite=TRUE)
message("Copied volcano plot")

# Table 5.3 — DMR count summary CSV
message("Table 5.3: DMR count summary...")
tbl <- do.call(rbind, lapply(CONTRASTS, function(ct) {
  d <- dmrs[[ct$name]]
  if (is.null(d)) return(NULL)
  data.frame(
    Contrast          = ct$label,
    Total_DMRs        = length(d),
    Hypomethylated    = sum(mcols(d)$regionType=="gain"),
    Hypermethylated   = sum(mcols(d)$regionType=="loss"),
    Median_size_bp    = round(median(width(d))),
    Median_CpG_count  = round(median(mcols(d)$cytosinesCount))
  )
}))
write.csv(tbl, file.path(OUT, "Table5.3_DMR_counts.csv"), row.names=FALSE)
print(tbl)
message("Saved Table5.3")

message("\nAll thesis figures saved to: ", OUT)
message("Files generated:")
for (f in list.files(OUT)) message("  ", f)
