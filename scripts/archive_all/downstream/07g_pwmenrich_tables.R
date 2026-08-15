#!/usr/bin/env Rscript
.libPaths('/data/home/bt25018/R/library')
library(ggplot2)

contrasts <- list(
  list(csv='results/tf_motif/pwmenrich_top_motifs_ASO_VPA_5k.csv',
       name='ASO_VPA_5k', title='ASO+VPA vs Scramble_CTRL (5k DMR set)'),
  list(csv='results/tf_motif/pwmenrich_top_motifs_VPA_5k.csv',
       name='VPA_5k', title='VPA vs Scramble_CTRL (5k DMR set)'),
  list(csv='results/tf_motif/pwmenrich_top_motifs_synergy_only_5k.csv',
       name='synergy_only_5k', title='Synergy-only / combination-emergent DMRs'),
  list(csv='results/tf_motif/pwmenrich_top_motifs_chr17_cacng.csv',
       name='CACNG_cluster', title='CACNG cluster locus-restricted (chr17:66.3-67.1 Mb)')
)

make_table_plot <- function(df, title) {
  df <- head(df, 20)
  df$target  <- gsub('Hsapiens-jolma2013-', '', df$target)
  df$p.value <- formatC(df$p.value, format='e', digits=2)
  df$raw.score <- sprintf('%.2f', df$raw.score)
  df$top.motif.prop <- paste0(round(df$top.motif.prop * 100), '%')

  # Build display table
  tbl <- data.frame(
    Rank            = as.character(seq_len(nrow(df))),
    Target          = df$target,
    PWM_Motif_ID    = df$id,
    Raw_score       = df$raw.score,
    P_value         = df$p.value,
    In_top_motifs   = df$top.motif.prop,
    stringsAsFactors = FALSE
  )

  # Column x positions (normalised 0-1)
  col_x     <- c(0.03, 0.10, 0.42, 0.62, 0.76, 0.92)
  col_names <- c('Rank', 'Target', 'PWM Motif ID', 'Raw score', 'P-value', 'In top motifs')
  n_rows    <- nrow(tbl)
  n_cols    <- length(col_x)

  # y positions: header at top, rows below
  header_y  <- n_rows + 1.5
  row_y     <- seq(n_rows, 1)   # row 1 at top of data

  # Build long data frame for geom_text
  cells <- data.frame()
  # Header row
  for (j in seq_len(n_cols)) {
    cells <- rbind(cells, data.frame(
      x=col_x[j], y=header_y,
      label=col_names[j], bold=TRUE, stringsAsFactors=FALSE))
  }
  # Data rows
  for (i in seq_len(n_rows)) {
    for (j in seq_len(n_cols)) {
      cells <- rbind(cells, data.frame(
        x=col_x[j], y=row_y[i],
        label=as.character(tbl[i, j]), bold=FALSE, stringsAsFactors=FALSE))
    }
  }

  # Alternating row background
  bg <- data.frame(
    ymin = seq(1, n_rows) - 0.5,
    ymax = seq(1, n_rows) + 0.5,
    fill = ifelse(seq(1, n_rows) %% 2 == 0, '#F7F7F7', 'white')
  )

  p <- ggplot() +
    # Alternating row shading
    geom_rect(data=bg,
              aes(xmin=0, xmax=1, ymin=ymin, ymax=ymax, fill=fill),
              colour=NA) +
    scale_fill_identity() +
    # Header background
    annotate('rect', xmin=0, xmax=1,
             ymin=header_y - 0.5, ymax=header_y + 0.5,
             fill='#2C3E50', colour=NA) +
    # Header separator line
    annotate('segment', x=0, xend=1,
             y=header_y - 0.5, yend=header_y - 0.5,
             colour='#2C3E50', linewidth=0.5) +
    # Cell text
    geom_text(data=cells[!cells$bold, ],
              aes(x=x, y=y, label=label),
              hjust=0, size=2.8, colour='grey15', family='mono') +
    # Header text
    geom_text(data=cells[cells$bold, ],
              aes(x=x, y=y, label=label),
              hjust=0, size=3.0, fontface='bold', colour='white') +
    # Title
    labs(title=title) +
    scale_x_continuous(limits=c(0, 1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0.4, header_y + 0.6), expand=c(0,0)) +
    theme_void() +
    theme(
      plot.title   = element_text(size=11, face='bold', margin=margin(b=8)),
      plot.margin  = margin(12, 12, 12, 12)
    )
  p
}

OUT_DIR <- 'results/tf_motif/plots'
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

for (cfg in contrasts) {
  cat('Plotting:', cfg$name, '\n')
  df <- read.csv(cfg$csv)
  p  <- make_table_plot(df, cfg$title)
  out <- file.path(OUT_DIR, paste0('pwmenrich_table_', cfg$name, '.pdf'))
  ggsave(out, p, width=11, height=8, device='pdf')
  cat('  Saved:', basename(out), '\n')
}
cat('Done\n')
