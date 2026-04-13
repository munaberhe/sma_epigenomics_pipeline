.libPaths("~/R/library")
# coverage_correlation.R
# WGBS QC Analysis: Coverage and Methylation Correlation
# SMA Epigenomics Pipeline — Muna Berhe, QMUL
# Run after all 12 methylation extractions are complete

library(data.table)
library(ggplot2)
library(corrplot)

# Config
COV_DIR <- "results/alignments/bs"
OUT_DIR <- "results/qc/methylation"
MIN_COV <- 5  # minimum coverage to include a CpG site

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Sample metadata
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

read_cov <- function(sample_name, min_cov = MIN_COV) {
  path <- file.path(COV_DIR,
    paste0(sample_name, "_bismark.deduplicated.bismark.cov.gz"))
  if (!file.exists(path)) {
    message("Missing: ", path)
    return(NULL)
  }
  message("Loading ", sample_name)
  dt <- fread(path, header = FALSE,
              col.names = c("chr", "start", "end",
                            "methylation", "count_m", "count_um"))
  dt[, coverage := count_m + count_um]
  dt[coverage >= min_cov]
}

cov_list <- lapply(SAMPLES$name, read_cov)
names(cov_list) <- SAMPLES$name

missing <- sapply(cov_list, is.null)
if (any(missing)) {
  message("Skipping missing samples: ",
          paste(names(missing)[missing], collapse = ", "))
  cov_list <- cov_list[!missing]
  SAMPLES  <- SAMPLES[SAMPLES$name %in% names(cov_list), ]
}

cov_summary <- do.call(rbind, lapply(names(cov_list), function(s) {
  dt <- cov_list[[s]]
  data.frame(
    sample        = s,
    group         = SAMPLES$group[SAMPLES$name == s],
    n_cpg_total   = nrow(dt),
    mean_coverage = round(mean(dt$coverage), 1),
    median_cov    = round(median(dt$coverage), 1),
    pct_5x        = round(mean(dt$coverage >= 5)  * 100, 1),
    pct_10x       = round(mean(dt$coverage >= 10) * 100, 1),
    global_meth   = round(mean(dt$methylation), 1),
    stringsAsFactors = FALSE
  )
}))

print(cov_summary)
write.csv(cov_summary, file.path(OUT_DIR, "coverage_summary.csv"),
          row.names = FALSE)

cov_plot_data <- do.call(rbind, lapply(names(cov_list), function(s) {
  dt <- cov_list[[s]]
  data.frame(
    sample   = s,
    group    = SAMPLES$group[SAMPLES$name == s],
    coverage = pmin(dt$coverage, 100)
  )
}))

p_cov <- ggplot(cov_plot_data, aes(x = coverage, colour = group)) +
  geom_density(alpha = 0.7, linewidth = 0.8) +
  facet_wrap(~sample, ncol = 4) +
  scale_colour_manual(values = GROUP_COLOURS) +
  geom_vline(xintercept = c(5, 10), linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  labs(title    = "CpG Coverage Distribution per Sample",
       subtitle = paste0("Dashed lines at 5x and 10x | Min coverage filter: ",
                         MIN_COV, "x"),
       x = "Coverage depth", y = "Density", colour = "Group") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7))

ggsave(file.path(OUT_DIR, "coverage_distribution.pdf"),
       p_cov, width = 12, height = 10)

cov_long <- reshape(
  cov_summary[, c("sample", "group", "pct_5x", "pct_10x")],
  varying   = c("pct_5x", "pct_10x"),
  v.names   = "pct",
  timevar   = "threshold",
  times     = c("5x coverage", "10x coverage"),
  direction = "long"
)

p_bar <- ggplot(cov_long,
                aes(x = sample, y = pct, fill = group, alpha = threshold)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = GROUP_COLOURS) +
  scale_alpha_manual(values = c("5x coverage" = 0.9, "10x coverage" = 0.5)) +
  geom_hline(yintercept = 80, linetype = "dashed",
             colour = "red", linewidth = 0.5) +
  labs(title    = "CpG Sites Meeting Coverage Thresholds",
       subtitle = "Red dashed line = 80% target",
       x = NULL, y = "% of CpG sites",
       fill = "Group", alpha = "Threshold") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "coverage_thresholds.pdf"),
       p_bar, width = 12, height = 6)

p_meth <- ggplot(cov_summary,
                 aes(x = group, y = global_meth,
                     colour = group, fill = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3) +
  scale_colour_manual(values = GROUP_COLOURS) +
  scale_fill_manual(values   = GROUP_COLOURS) +
  geom_hline(yintercept = c(70, 80), linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  labs(title    = "Global CpG Methylation per Sample",
       subtitle = "Expected: 70-80% for human somatic cells",
       x = "Group", y = "Mean CpG methylation (%)",
       colour = "Group", fill = "Group") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUT_DIR, "global_methylation.pdf"),
       p_meth, width = 8, height = 6)

make_key <- function(dt) paste0(dt$chr, ":", dt$start)
all_keys  <- Reduce(intersect, lapply(cov_list, make_key))
message("Common CpG sites: ", format(length(all_keys), big.mark = ","))

if (length(all_keys) < 1000) {
  warning("Too few common sites — check coverage filters")
} else {
  meth_matrix <- do.call(cbind, lapply(names(cov_list), function(s) {
    dt  <- cov_list[[s]]
    key <- make_key(dt)
    dt[match(all_keys, key), methylation]
  }))
  colnames(meth_matrix) <- names(cov_list)

  cor_matrix <- cor(meth_matrix, use = "complete.obs", method = "pearson")
  write.csv(round(cor_matrix, 4),
            file.path(OUT_DIR, "methylation_correlation.csv"))

  label_colours <- GROUP_COLOURS[
    SAMPLES$group[match(colnames(cor_matrix), SAMPLES$name)]
  ]

  pdf(file.path(OUT_DIR, "methylation_correlation_heatmap.pdf"),
      width = 10, height = 9)
  corrplot(cor_matrix,
           method      = "color",
           type        = "upper",
           order       = "hclust",
           col         = colorRampPalette(c("#065A82", "white", "#F59E0B"))(200),
           addCoef.col = "black",
           number.cex  = 0.65,
           tl.col      = label_colours,
           tl.cex      = 0.8,
           tl.srt      = 45,
           cl.lim      = c(0.85, 1),
           title       = "Pairwise Methylation Correlation (Pearson)",
           mar         = c(0, 0, 2, 0))
  dev.off()

  low_cor <- which(cor_matrix < 0.95 & upper.tri(cor_matrix), arr.ind = TRUE)
  if (nrow(low_cor) > 0) {
    message("Low correlations (<0.95) detected:")
    for (i in seq_len(nrow(low_cor))) {
      r <- low_cor[i, 1]
      c <- low_cor[i, 2]
      message("  ", rownames(cor_matrix)[r], " vs ",
              colnames(cor_matrix)[c], ": ",
              round(cor_matrix[r, c], 4))
    }
  } else {
    message("All pairwise correlations >= 0.95")
  }
}

message("Done. Outputs in: ", OUT_DIR)
