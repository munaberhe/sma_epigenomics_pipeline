
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages(library(PWMEnrich))
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/tf_motif"

BLACKLIST <- "CENPB|ZNF274|ZNF93|^UW\\.Motif"

# Load pre-computed CSVs and plot directly - no re-running motifEnrichment
CONTRASTS <- list(
  list(name="VPA",          csv="pwmenrich_top_motifs_VPA.csv"),
  list(name="ASO_VPA",      csv="pwmenrich_top_motifs_ASO_VPA.csv"),
  list(name="synergy_only", csv="pwmenrich_top_motifs_synergy_only.csv")
)

for (ct in CONTRASTS) {
  path <- file.path(OUT_DIR, ct$csv)
  if (!file.exists(path)) { message("missing: ", path); next }
  message("=== ", ct$name, " ===")
  
  df <- read.csv(path, stringsAsFactors=FALSE)
  df <- df[!grepl(BLACKLIST, df$target, ignore.case=TRUE, perl=TRUE),]
  df <- df[order(df$p.value),]
  df <- df[!duplicated(df$target),]  # remove duplicate TF names
  
  # top 20 report - use ggplot bar chart since we have CSV not RDS
  library(ggplot2)
  top20 <- head(df, 20)
  top20$target <- factor(top20$target, levels=rev(top20$target))
  top20$neg_logp <- -log10(top20$p.value + 1e-300)
  
  p <- ggplot(top20, aes(x=target, y=neg_logp)) +
    geom_col(fill="#1F3A5F") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed",
               colour="#C0392B", linewidth=0.5) +
    coord_flip() +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold")) +
    labs(title=paste0("Top 20 enriched TF motifs: ", ct$name),
         subtitle="PWMEnrich MotifDb background | dashed = p=0.05 threshold",
         x=NULL, y=expression(-log[10](p)))
  
  # top 20
  fname <- file.path(OUT_DIR, paste0("pwmenrich_report_", ct$name, ".pdf"))
  ggsave(fname, p, width=10, height=8, device=cairo_pdf)
  message("saved: ", basename(fname))
  
  # top 10
  top10 <- head(df, 10)
  top10$target <- factor(top10$target, levels=rev(top10$target))
  top10$neg_logp <- -log10(top10$p.value + 1e-300)
  p10 <- ggplot(top10, aes(x=target, y=neg_logp)) +
    geom_col(fill="#1F3A5F") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed",
               colour="#C0392B", linewidth=0.5) +
    coord_flip() +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold")) +
    labs(title=paste0("Top 10 enriched TF motifs: ", ct$name),
         subtitle="PWMEnrich MotifDb background",
         x=NULL, y=expression(-log[10](p)))
  
  fname2 <- file.path(OUT_DIR, paste0("pwmenrich_top10_", ct$name, ".pdf"))
  ggsave(fname2, p10, width=9, height=6, device=cairo_pdf)
  message("saved: ", basename(fname2))
}
message("Done.")
