#!/usr/bin/env Rscript
.libPaths('/data/home/bt25018/R/library')
suppressPackageStartupMessages({ library(PWMEnrich); library(ggplot2) })

OUT_DIR <- 'results/tf_motif/plots'
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

extract_top <- function(report, n=15) {
  df <- report@d
  pval_col <- grep('p.val|pval|p_val', names(df), value=TRUE, ignore.case=TRUE)[1]
  df$pval <- df[[pval_col]]
  df <- head(df, n)
  df$label <- gsub('Hsapiens-jolma2013-', '', df$target)
  df$label <- gsub('Hsapiens-', '', df$label)
  # Deduplicate: any repeated label gets motif ID appended
  seen <- c()
  for (i in seq_len(nrow(df))) {
    lbl <- df$label[i]
    if (lbl %in% seen) {
      df$label[i] <- paste0(lbl, '_', i)
    }
    seen <- c(seen, df$label[i])
  }
  df$neglog10p <- -log10(df$pval)
  df$label <- factor(df$label, levels=rev(df$label))
  df
}

REPORTS <- list(
  list(rds='results/tf_motif/pwmenrich_report.rds',
       name='ASO_specific_151dmrs', color='#02C39A',
       title='ASO-specific DMRs (n=151)'),
  list(rds='results/tf_motif/pwmenrich_report_ASO_VPA_5k.rds',
       name='ASO_VPA_vs_Scramble_CTRL_5k', color='#F59E0B',
       title='ASO+VPA vs Scramble_CTRL (5k DMR set)'),
  list(rds='results/tf_motif/pwmenrich_report_VPA_5k.rds',
       name='VPA_vs_Scramble_CTRL_5k', color='#1C7293',
       title='VPA vs Scramble_CTRL (5k DMR set)'),
  list(rds='results/tf_motif/pwmenrich_report_synergy_only_5k.rds',
       name='synergy_only_5k', color='#D4A017',
       title='Synergy-only / combination-emergent DMRs')
)

for (cfg in REPORTS) {
  cat('Plotting:', cfg$name, '\n')
  if (!file.exists(cfg$rds)) { cat('  SKIP: file not found\n'); next }
  report <- readRDS(cfg$rds)
  df <- tryCatch(
    extract_top(report),
    error = function(e) { cat('  ERROR:', e$message, '\n'); NULL }
  )
  if (is.null(df)) next
  p <- ggplot(df, aes(x=neglog10p, y=label)) +
    geom_col(fill=cfg$color, width=0.7) +
    geom_text(aes(label=sprintf('p=%.1e', pval)),
              hjust=-0.05, size=3.0, colour='grey30') +
    scale_x_continuous(expand=expansion(mult=c(0, 0.25)),
                       name='-log10(p-value)') +
    labs(title=cfg$title, y=NULL) +
    theme_bw(base_size=11) +
    theme(
      plot.title         = element_text(size=11, face='bold'),
      axis.text.y        = element_text(size=9, face='bold'),
      axis.text.x        = element_text(size=8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(10, 20, 10, 10)
    )
  outfile <- file.path(OUT_DIR, paste0(cfg$name, '_barplot.pdf'))
  ggsave(outfile, p, width=7, height=5.5, device='pdf')
  cat('  Saved:', basename(outfile), '\n')
}

# Combined 4-panel
cat('\nBuilding combined plot...\n')
panels <- list()
for (cfg in REPORTS) {
  if (!file.exists(cfg$rds)) next
  report <- readRDS(cfg$rds)
  df <- tryCatch(extract_top(report, n=10), error=function(e) NULL)
  if (is.null(df)) next
  p <- ggplot(df, aes(x=neglog10p, y=label)) +
    geom_col(fill=cfg$color, width=0.7) +
    scale_x_continuous(expand=expansion(mult=c(0, 0.2))) +
    labs(title=cfg$title, x='-log10(p)', y=NULL) +
    theme_bw(base_size=9) +
    theme(
      plot.title         = element_text(size=8, face='bold'),
      axis.text.y        = element_text(size=7.5),
      axis.text.x        = element_text(size=7),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(6, 12, 6, 6)
    )
  panels[[length(panels)+1]] <- p
}

outfile <- file.path(OUT_DIR, 'all_contrasts_combined.pdf')
if (requireNamespace('patchwork', quietly=TRUE)) {
  library(patchwork)
  combined <- Reduce('+', panels) + plot_layout(ncol=2)
  ggsave(outfile, combined, width=12, height=10, device='pdf')
} else {
  library(gridExtra)
  pdf(outfile, width=12, height=10)
  do.call(grid.arrange, c(panels, ncol=2))
  dev.off()
}
cat('  Saved: all_contrasts_combined.pdf\n')
cat('\nDone. Outputs in:', OUT_DIR, '\n')
