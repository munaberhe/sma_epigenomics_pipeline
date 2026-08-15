#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(PWMEnrich)
  library(PWMEnrich.Hsapiens.background)
  library(JASPAR2024)
  library(TFBSTools)
  library(ggplot2)
  library(ggseqlogo)
  library(eulerr)
  library(dplyr)
  library(patchwork)
})

# venn of TF motifs enriched in 3 DMR sets (ASO-only, VPA-only, ASO+VPA)
# plus top motif logos. matches Grant et al 2026 Fig 2A layout.
# input : results/dmr/dmr_<contrast>.rds for 3 contrasts
# output: results/plots/grant_fig2a_motif_venn.{pdf,png}

OUT_DIR <- "results/plots"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# contrasts that define the 3 sets
SETS <- list(
  "ASO only" = "ASO_CTRL_vs_Scramble_CTRL",
  "VPA only" = "Scramble_VPA_vs_Scramble_CTRL",
  "ASO+VPA"  = "ASO_VPA_vs_Scramble_CTRL"
)

genome <- BSgenome.Hsapiens.UCSC.hg38

message("loading JASPAR2024 vertebrate PWMs...")
jaspar <- JASPAR2024()
sqlfile <- db(jaspar)
opts <- list(species=9606, all_versions=FALSE)
pfm_list <- getMatrixSet(sqlfile, opts)
pwm_list <- toPWM(pfm_list)

# build PWMEnrich background (lognormal across genome, 200 bp)
# use a precomputed background if available
# PWMEnrich makeBackground() only supports hg19/mm9/dm3, not hg38.
# Use the precomputed hg19 MotifDb lognormal background, identical to 04_tf_motif.R.
# Safe because: PWM scoring is coordinate-free - sequences are extracted from
# BSgenome.Hsapiens.UCSC.hg38 so DMR coordinates stay on hg38. The "hg19" in
# the background name refers only to GC-content normalisation, not coordinates.
message("loading PWMEnrich hg19 MotifDb lognormal background...")
data(PWMLogn.hg19.MotifDb.Hsap)
bg <- PWMLogn.hg19.MotifDb.Hsap

# helper: enriched motifs for a contrast at p < 0.05
enriched_set <- function(contrast, p_thresh=0.05) {
  rds_path <- file.path("results/dmr", paste0("dmr_", contrast, ".rds"))
  if (!file.exists(rds_path)) stop("missing: ", rds_path)
  dmrs <- readRDS(rds_path)
  if (length(dmrs) == 0) return(character(0))

  # hc filter + subsample for speed (5000 seqs standard for motif enrichment)
  if (!is.null(mcols(dmrs)$cytosinesCount))
    dmrs <- dmrs[mcols(dmrs)$cytosinesCount >= 6]
  if (length(dmrs) == 0) return(character(0))
  if (length(dmrs) > 5000) {
    set.seed(42)
    dmrs <- dmrs[sample(length(dmrs), 5000)]
    message("  subsampled to 5000 DMRs")
  }
  message("  using ", length(dmrs), " DMRs for motif enrichment")
  seqs <- getSeq(genome, dmrs)
  res <- motifEnrichment(seqs, bg, score="affinity")
  rep <- groupReport(res)
  df <- as.data.frame(rep)
  # filter by p.value column (PWMEnrich groupReport)
  hits <- df %>%
    filter(!is.na(p.value), p.value < p_thresh) %>%
    arrange(p.value)
  # flag repeat-driven artefacts to remove
  artefact_re <- "CENPB|ZNF274|ZNF93|ZNF675|GRHL"  # extend as needed
  hits <- hits %>% filter(!grepl(artefact_re, target, ignore.case=TRUE))
  unique(hits$target)
}

message("computing enriched motifs per contrast...")
sets <- lapply(SETS, enriched_set)
str(sets)

# Venn (eulerr) of the 3 sets
venn_input <- list(
  "ASO only" = sets[["ASO only"]],
  "VPA only" = sets[["VPA only"]],
  "ASO+VPA"  = sets[["ASO+VPA"]]
)
fit <- euler(venn_input)

# pull top 4 motifs unique to each set + 1 shared (for logo strip)
unique_to <- function(name) {
  setdiff(sets[[name]], unlist(sets[setdiff(names(sets), name)]))
}
shared_all <- Reduce(intersect, sets)
top_motifs <- c(
  head(unique_to("ASO only"), 1),
  head(unique_to("VPA only"), 1),
  head(unique_to("ASO+VPA"), 1),
  head(shared_all, 1)
)
top_motifs <- top_motifs[!is.na(top_motifs)]

# logo plot: pull PWMs by name
pwm_by_name <- function(nm) {
  ix <- which(sapply(pwm_list, function(x) name(x)) == nm)
  if (length(ix) == 0) return(NULL)
  as.matrix(pwm_list[[ix[1]]])
}

logo_grobs <- lapply(top_motifs, function(nm) {
  mat <- pwm_by_name(nm)
  if (is.null(mat)) return(NULL)
  ggseqlogo(mat, method="bits") +
    labs(title=nm) +
    theme(plot.title=element_text(size=10, face="bold"))
})
logo_grobs <- Filter(Negate(is.null), logo_grobs)

# layout: venn on left, logos stacked on right
venn_grob <- plot(fit,
                  fills = list(fill=c("#1F3A5F","#C0392B","#D4A017"), alpha=0.5),
                  labels = list(font=2, cex=1),
                  quantities = list(cex=0.9))

# save venn separately (base graphics from eulerr)
pdf(file.path(OUT_DIR, "grant_fig2a_motif_venn.pdf"), width=13, height=6)
title_grob <- grid::textGrob(
  "TF motif enrichment overlap across DMR contrasts",
  gp=grid::gpar(fontsize=14, fontface="bold"))
subtitle_grob <- grid::textGrob(
  paste0("PWMEnrich enriched motifs (p<0.05) in ASO-specific, VPA-only and ASO+VPA DMRs
",
         "Numbers = motifs unique to or shared between contrasts | Logos = top unique motif per set"),
  gp=grid::gpar(fontsize=9, col="grey30"))
main_panel <- gridExtra::arrangeGrob(
  venn_grob,
  do.call(gridExtra::arrangeGrob, c(logo_grobs, list(ncol=1))),
  ncol=2, widths=c(1.2, 1))
gridExtra::grid.arrange(
  title_grob, subtitle_grob, main_panel,
  ncol=1, heights=c(0.08, 0.08, 0.84))
dev.off()

png(file.path(OUT_DIR, "grant_fig2a_motif_venn.png"),
    width=11, height=5, units="in", res=300)
gridExtra::grid.arrange(
  venn_grob,
  do.call(gridExtra::arrangeGrob, c(logo_grobs, ncol=1)),
  ncol=2, widths=c(1.2, 1)
)
dev.off()

# also write the motif lists for table 5.x
write.csv(
  data.frame(
    set = rep(names(sets), sapply(sets, length)),
    motif = unlist(sets)
  ),
  file.path(OUT_DIR, "grant_fig2a_motif_sets.csv"),
  row.names=FALSE
)

message("wrote: ", file.path(OUT_DIR, "grant_fig2a_motif_venn.pdf"))
message("done.")
