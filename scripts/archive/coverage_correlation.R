.libPaths("~/R/library")
# coverage_correlation.R
# WGBS QC Analysis: Coverage and Methylation Correlation
# SMA Epigenomics Pipeline — Muna Berhe, QMUL
# CpG-only version: shell-level join filtering to minimise memory usage

library(data.table)
library(ggplot2)
library(corrplot)

COV_DIR  <- "results/alignments/bs"
OUT_DIR  <- "results/qc/methylation"
MIN_COV  <- 5
CPG_POS  <- file.path(COV_DIR, "cpg_positions_hg38.txt")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- data.frame(
  name  = c("Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
            "Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3",
            "ASO_CTRL_1",      "ASO_CTRL_2",      "ASO_CTRL_3",
            "ASO_VPA_1",       "ASO_VPA_2",       "ASO_VPA_3"),
  group = c(rep("Scramble_CTRL", 3), rep("Scramble_VPA", 3),
            rep("ASO_CTRL", 3),      rep("ASO_VPA", 3)),
  stringsAsFactors = FALSE
)

GROUP_COLOURS <- c(
  "Scramble_CTRL" = "#065A82",
  "Scramble_VPA"  = "#1C7293",
  "ASO_CTRL"      = "#02C39A",
  "ASO_VPA"       = "#F59E0B"
)

# Helper: shell command to filter cov.gz to CpG-only sites via join
# For full genome (chrom=NULL) or single chromosome
make_cmd <- function(path, chrom = NULL) {
  if (is.null(chrom)) {
    pos_cmd <- paste0("cat ", CPG_POS)
  } else {
    pos_cmd <- paste0("grep '^", chrom, "\t' ", CPG_POS)
  }
  paste0(
    "join ",
    "<(", pos_cmd, " | awk '{print $1\"_\"$2}' | sort) ",
    "<(zcat ", path, " | awk '{print $1\"_\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6}' | sort -k1,1) ",
    "| awk '{split($1,a,\"_\"); print a[1]\"\\t\"a[2]\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5}'"
  )
}

# Per-sample coverage statistics (CpG-only)
message("Computing per-sample coverage statistics (CpG-only)...")
cov_summary_list <- list()

for (s in SAMPLES$name) {
  path <- file.path(COV_DIR, paste0(s, "_bismark.deduplicated.bismark.cov.gz"))
  if (!file.exists(path)) { message("Missing: ", path); next }
  message("Processing ", s)
  dt <- fread(cmd = make_cmd(path),
              header = FALSE,
              col.names = c("chr", "start", "end", "methylation", "count_m", "count_um"))
  dt[, coverage := count_m + count_um]
  dt <- dt[coverage >= MIN_COV]
  cov_summary_list[[s]] <- data.frame(
    sample        = s,
    group         = SAMPLES$group[SAMPLES$name == s],
    n_cpg_total   = nrow(dt),
    mean_coverage = round(mean(dt$coverage), 1),
    median_cov    = round(median(dt$coverage), 1),
    pct_5x        = round(mean(dt$coverage >= 5)  * 100, 1),
    pct_10x       = round(mean(dt$coverage >= 10) * 100, 1),
    global_meth   = round(sum(dt$count_m) /
                          sum(dt$count_m + dt$count_um) * 100, 1),
    stringsAsFactors = FALSE
  )
  rm(dt); gc()
}

cov_summary <- do.call(rbind, cov_summary_list)
print(cov_summary)
write.csv(cov_summary, file.path(OUT_DIR, "coverage_summary.csv"), row.names = FALSE)
message("Saved: coverage_summary.csv")

# Coverage distribution plots
message("Plotting coverage distributions...")
cov_plot_list <- list()
for (s in SAMPLES$name) {
  path <- file.path(COV_DIR, paste0(s, "_bismark.deduplicated.bismark.cov.gz"))
  if (!file.exists(path)) next
  message("Loading for plot: ", s)
  dt <- fread(cmd = make_cmd(path),
              header = FALSE,
              col.names = c("chr", "start", "end", "methylation", "count_m", "count_um"))
  dt[, coverage := count_m + count_um]
  dt <- dt[coverage >= MIN_COV]
  n   <- min(500000, nrow(dt))
  idx <- sample(nrow(dt), n)
  cov_plot_list[[s]] <- data.frame(
    sample   = s,
    group    = SAMPLES$group[SAMPLES$name == s],
    coverage = pmin(dt$coverage[idx], 100)
  )
  rm(dt); gc()
}
cov_plot_data <- do.call(rbind, cov_plot_list)
rm(cov_plot_list); gc()

p_cov <- ggplot(cov_plot_data, aes(x = coverage, colour = group)) +
  geom_density(alpha = 0.7, linewidth = 0.8) +
  facet_wrap(~sample, ncol = 4) +
  scale_colour_manual(values = GROUP_COLOURS) +
  geom_vline(xintercept = c(5, 10), linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  labs(title    = "CpG Coverage Distribution per Sample",
       subtitle = paste0("Dashed lines at 5x and 10x | Min coverage: ", MIN_COV, "x"),
       x = "Coverage depth", y = "Density", colour = "Group") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))
ggsave(file.path(OUT_DIR, "coverage_distribution.pdf"), p_cov, width = 12, height = 10)
message("Saved: coverage_distribution.pdf")
rm(cov_plot_data, p_cov); gc()

# Coverage threshold bar chart
cov_long <- reshape(
  cov_summary[, c("sample", "group", "pct_5x", "pct_10x")],
  varying   = c("pct_5x", "pct_10x"),
  v.names   = "pct",
  timevar   = "threshold",
  times     = c("5x coverage", "10x coverage"),
  direction = "long"
)
p_bar <- ggplot(cov_long, aes(x = sample, y = pct, fill = group, alpha = threshold)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = GROUP_COLOURS) +
  scale_alpha_manual(values = c("5x coverage" = 0.9, "10x coverage" = 0.5)) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "red", linewidth = 0.5) +
  labs(title    = "CpG Sites Meeting Coverage Thresholds",
       subtitle = "Red dashed line = 80% target",
       x = NULL, y = "% of CpG sites", fill = "Group", alpha = "Threshold") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggsave(file.path(OUT_DIR, "coverage_thresholds.pdf"), p_bar, width = 12, height = 6)
message("Saved: coverage_thresholds.pdf")
rm(p_bar); gc()

# Global methylation plot
p_meth <- ggplot(cov_summary, aes(x = group, y = global_meth,
                                   colour = group, fill = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3) +
  scale_colour_manual(values = GROUP_COLOURS) +
  scale_fill_manual(values = GROUP_COLOURS) +
  geom_hline(yintercept = c(70, 80), linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  labs(title    = "Global CpG Methylation per Sample",
       subtitle = "Expected: 70-80% for human somatic cells",
       x = "Group", y = "Mean CpG methylation (%)",
       colour = "Group", fill = "Group") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "global_methylation.pdf"), p_meth, width = 8, height = 6)
message("Saved: global_methylation.pdf")
rm(p_meth); gc()

# Methylation correlation (chr1, CpG-only)
message("Computing methylation correlation using chr1 CpG subset...")
meth_list <- list()
for (s in SAMPLES$name) {
  path <- file.path(COV_DIR, paste0(s, "_bismark.deduplicated.bismark.cov.gz"))
  if (!file.exists(path)) next
  message("Loading chr1 for correlation: ", s)
  dt <- fread(cmd = make_cmd(path, chrom = "chr1"),
              header = FALSE,
              col.names = c("chr", "start", "end", "methylation", "count_m", "count_um"))
  dt[, coverage := count_m + count_um]
  dt <- dt[coverage >= MIN_COV]
  dt[, key := paste0(chr, ":", start)]
  meth_list[[s]] <- dt[, .(key, methylation)]
  rm(dt); gc()
}

all_keys <- Reduce(intersect, lapply(meth_list, function(x) x$key))
message("Common CpG sites on chr1: ", format(length(all_keys), big.mark = ","))

meth_matrix <- do.call(cbind, lapply(names(meth_list), function(s) {
  m <- meth_list[[s]]
  m$methylation[match(all_keys, m$key)]
}))
colnames(meth_matrix) <- names(meth_list)
rm(meth_list); gc()

cor_matrix <- cor(meth_matrix, use = "complete.obs", method = "pearson")
write.csv(round(cor_matrix, 4), file.path(OUT_DIR, "methylation_correlation.csv"))
message("Saved: methylation_correlation.csv")

label_colours <- GROUP_COLOURS[SAMPLES$group[match(colnames(cor_matrix), SAMPLES$name)]]

pdf(file.path(OUT_DIR, "methylation_correlation_heatmap.pdf"), width = 10, height = 9)
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
         col          = colorRampPalette(c("#065A82", "white", "#F59E0B"))(200),
         addCoef.col  = "black", number.cex = 0.65,
         tl.col       = label_colours, tl.cex = 0.8, tl.srt = 45,
         cl.lim       = c(0.85, 1),
         title        = "Pairwise Methylation Correlation (Pearson) — chr1 CpG",
         mar          = c(0, 0, 2, 0))
dev.off()
message("Saved: methylation_correlation_heatmap.pdf")

low_cor <- which(cor_matrix < 0.95 & upper.tri(cor_matrix), arr.ind = TRUE)
if (nrow(low_cor) > 0) {
  message("Low correlations (<0.95) detected:")
  for (i in seq_len(nrow(low_cor))) {
    r <- low_cor[i, 1]; c <- low_cor[i, 2]
    message("  ", rownames(cor_matrix)[r], " vs ", colnames(cor_matrix)[c],
            ": ", round(cor_matrix[r, c], 4))
  }
} else {
  message("All pairwise correlations >= 0.95")
}

message("Done. Outputs in: ", OUT_DIR)
