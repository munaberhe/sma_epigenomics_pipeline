#!/usr/bin/env Rscript
# 17_pairwise_context_scan.R
# Identifies context-dependent and pairwise synergy candidate genes from
# the four pairwise DMR contrasts.
#
# Outputs:
#   results/pairwise_context_scan/ASO_context_dependent_scored.csv
#   results/pairwise_context_scan/VPA_context_dependent_scored.csv
#   results/pairwise_context_scan/all_pairwise_candidates.csv
#   results/pairwise_context_scan/pairwise_synergy_candidates_provenance.csv
#   results/pairwise_context_scan/appendix_table_F1.csv
#
# Context-dependent definition:
#   ASO context-dependent: DMR called in ASO_VPA_vs_Scramble_VPA (ASO on VPA
#     background) but NOT in ASO_CTRL_vs_Scramble_CTRL (ASO alone)
#   VPA context-dependent: DMR called in ASO_VPA_vs_ASO_CTRL (VPA on ASO
#     background) but NOT in Scramble_VPA_vs_Scramble_CTRL (VPA alone)
#
# Pairwise synergy candidates:
#   Genes appearing in BOTH the ASO context-dependent and VPA context-dependent
#   lists, indicating their methylation is altered by the combination in a
#   context-dependent way under both treatment backgrounds.
#
# Relevance scoring (0-10):
#   +3 promoter or UTR annotation
#   +2 known SMA/motor neuron gene (manual flag)
#   +2 p < 0.001
#   +2 meth_diff > 0.30
#   +1 CpG island overlap

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/pairwise_context_scan"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# load DMR results
message("Loading DMR results...")
dmr <- list()
contrasts <- c(
  "ASO_CTRL_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_VPA_vs_Scramble_VPA",
  "ASO_VPA_vs_ASO_CTRL"
)
for (ct in contrasts) {
  rds <- paste0("results/dmr/dmr_", ct, ".rds")
  if (file.exists(rds)) {
    dmr[[ct]] <- readRDS(rds)
    message("  ", ct, ": ", length(dmr[[ct]]))
  }
}

# load annotation for the two context contrasts
load_annotated <- function(contrast) {
  csv <- paste0("results/dmr_annotation/", contrast, "_annotated.csv")
  if (!file.exists(csv)) {
    message("  Annotation CSV not found: ", csv)
    return(NULL)
  }
  df <- read.csv(csv, stringsAsFactors=FALSE)
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  df
}

message("Loading annotated DMR tables...")
aso_context <- load_annotated("ASO_VPA_vs_Scramble_VPA")
vpa_context  <- load_annotated("ASO_VPA_vs_ASO_CTRL")
aso_alone    <- load_annotated("ASO_CTRL_vs_Scramble_CTRL")
vpa_alone    <- load_annotated("Scramble_VPA_vs_Scramble_CTRL")

# genes called in ASO alone / VPA alone - used for exclusion
aso_alone_genes <- if (!is.null(aso_alone)) unique(aso_alone$SYMBOL) else character(0)
vpa_alone_genes <- if (!is.null(vpa_alone)) unique(vpa_alone$SYMBOL) else character(0)

# SMA/motor neuron relevant genes (manually curated)
sma_relevant <- c(
  "SMN1","SMN2","NAIP","NCALD","PLS3","NRXN1","NRXN3",
  "CAMK2A","EPHB1","SHANK2","KIF21B","DSCAM","SYT1",
  "ROBO1","ROBO2","SLIT1","SLIT2","NRG1","ERBB4",
  "KDM1A","LSD1","HDAC1","HDAC2","DNMT1","DNMT3A","DNMT3B",
  "PAX5","IRF8","USP7","USP27X","ZDHHC22","AFAP1L1",
  "UBE2D4","TNS4","RELL2","DDIT4L","MRPS2","GNG14",
  "RNA5S13","KIAA1656","TCEAL4"
)

# relevance scoring function
score_gene <- function(df_row) {
  score <- 0
  ann <- df_row$annotation
  if (!is.null(ann) && !is.na(ann)) {
    if (grepl("Promoter|UTR", ann)) score <- score + 3
  }
  sym <- df_row$SYMBOL
  if (!is.null(sym) && !is.na(sym) && sym %in% sma_relevant) score <- score + 2
  pv <- as.numeric(df_row$pValue)
  if (!is.na(pv) && pv < 0.001) score <- score + 2
  md <- abs(as.numeric(df_row$meth_diff))
  if (!is.na(md) && md > 0.30) score <- score + 2
  # CpG island overlap placeholder - add +1 if annotation mentions CpG
  if (!is.null(ann) && !is.na(ann) && grepl("CpG", ann)) score <- score + 1
  score
}

# build scored table for a context-dependent set
build_scored <- function(df, exclude_genes, group_label) {
  if (is.null(df)) return(NULL)
  # keep genes NOT called in the single-treatment contrast
  df_ctx <- df[!df$SYMBOL %in% exclude_genes, ]
  if (nrow(df_ctx) == 0) return(NULL)
  # keep best DMR per gene (highest meth_diff)
  df_ctx$meth_diff <- abs(as.numeric(df_ctx$proportion2) -
                          as.numeric(df_ctx$proportion1))
  df_ctx <- df_ctx[order(-df_ctx$meth_diff), ]
  df_best <- df_ctx[!duplicated(df_ctx$SYMBOL), ]
  df_best$score <- apply(df_best, 1, score_gene)
  df_best$sma_relevant <- df_best$SYMBOL %in% sma_relevant
  df_best$group <- group_label
  df_best$effect_rank <- rank(-df_best$meth_diff) / nrow(df_best)
  df_best[order(-df_best$score, -df_best$meth_diff), ]
}

message("Building context-dependent candidate lists...")
aso_scored <- build_scored(aso_context, aso_alone_genes, "ASO_context_dependent")
vpa_scored  <- build_scored(vpa_context,  vpa_alone_genes,  "VPA_context_dependent")

if (!is.null(aso_scored)) {
  write.csv(aso_scored, file.path(OUT, "ASO_context_dependent_scored.csv"),
            row.names=FALSE)
  message("  ASO context-dependent: ", nrow(aso_scored), " genes")
}
if (!is.null(vpa_scored)) {
  write.csv(vpa_scored, file.path(OUT, "VPA_context_dependent_scored.csv"),
            row.names=FALSE)
  message("  VPA context-dependent: ", nrow(vpa_scored), " genes")
}

# pairwise synergy candidates - genes in both lists
message("Identifying pairwise synergy candidates...")
if (!is.null(aso_scored) && !is.null(vpa_scored)) {
  both_genes <- intersect(aso_scored$SYMBOL, vpa_scored$SYMBOL)
  message("  Genes in both lists: ", length(both_genes))

  # build provenance table
  prov <- do.call(rbind, lapply(both_genes, function(sym) {
    a <- aso_scored[aso_scored$SYMBOL == sym, ][1, ]
    v <- vpa_scored[vpa_scored$SYMBOL == sym, ][1, ]
    combined_score <- as.integer(a$score) + as.integer(v$score)
    data.frame(
      SYMBOL         = sym,
      GENENAME       = a$GENENAME,
      chr            = a$seqnames,
      dmr_start      = a$start,
      dmr_end        = a$end,
      ASO_annotation = a$annotation,
      VPA_annotation = v$annotation,
      ASO_score      = a$score,
      VPA_score      = v$score,
      combined_score = combined_score,
      ASO_meth_diff  = round(a$meth_diff, 3),
      VPA_meth_diff  = round(v$meth_diff, 3),
      geneStart      = a$geneStart,
      geneEnd        = a$geneEnd,
      strand         = a$geneStrand,
      selected       = ifelse(combined_score >= 4, "YES", "NO"),
      selection_reason = ifelse(combined_score >= 4,
        paste0("combined_score=", combined_score),
        paste0("combined_score=", combined_score, ", below threshold")),
      stringsAsFactors = FALSE
    )
  }))
  prov <- prov[order(-prov$combined_score), ]
  write.csv(prov, file.path(OUT, "pairwise_synergy_candidates_provenance.csv"),
            row.names=FALSE)
  message("  Written provenance CSV")

  # all pairwise candidates combined
  all_genes <- rbind(
    aso_scored[, c("SYMBOL","GENENAME","group","seqnames","start","end",
                   "annotation","distanceToTSS","score","meth_diff",
                   "regionType","sma_relevant","effect_rank")],
    vpa_scored[!vpa_scored$SYMBOL %in% aso_scored$SYMBOL,
               c("SYMBOL","GENENAME","group","seqnames","start","end",
                 "annotation","distanceToTSS","score","meth_diff",
                 "regionType","sma_relevant","effect_rank")]
  )
  # mark synergy genes
  all_genes$group[all_genes$SYMBOL %in% both_genes] <- "ASO+VPA_context_dependent"
  write.csv(all_genes, file.path(OUT, "all_pairwise_candidates.csv"),
            row.names=FALSE)
  message("  Written all_pairwise_candidates.csv: ", nrow(all_genes), " genes")

  # Appendix Table F.1 - clean formatted version for thesis
  appendix <- prov[prov$selected == "YES", ]
  appendix_ctx <- rbind(
    aso_scored[!aso_scored$SYMBOL %in% both_genes &
               aso_scored$score >= 4,
               c("SYMBOL","GENENAME","seqnames","start","end",
                 "annotation","score","meth_diff","group")],
    vpa_scored[!vpa_scored$SYMBOL %in% both_genes &
               vpa_scored$score >= 4,
               c("SYMBOL","GENENAME","seqnames","start","end",
                 "annotation","score","meth_diff","group")]
  )
  # format appendix table
  f1 <- data.frame(
    Category = c(
      rep("Pairwise synergy (ASO+VPA context-dependent)", nrow(appendix)),
      rep("ASO context-restricted", sum(!appendix_ctx$group %in%
          c("VPA_context_dependent","ASO+VPA_context_dependent"))),
      rep("VPA context-restricted", sum(appendix_ctx$group ==
          "VPA_context_dependent"))
    ),
    Symbol      = c(appendix$SYMBOL, appendix_ctx$SYMBOL),
    Gene_name   = c(appendix$GENENAME, appendix_ctx$GENENAME),
    Chr         = c(appendix$chr, appendix_ctx$seqnames),
    DMR_start   = c(appendix$dmr_start, appendix_ctx$start),
    DMR_end     = c(appendix$dmr_end, appendix_ctx$end),
    Annotation  = c(appendix$ASO_annotation, appendix_ctx$annotation),
    Score       = c(appendix$combined_score, appendix_ctx$score),
    Meth_diff   = c(appendix$ASO_meth_diff, appendix_ctx$meth_diff),
    stringsAsFactors = FALSE
  )
  write.csv(f1, file.path(OUT, "appendix_table_F1.csv"), row.names=FALSE)
  message("  Written appendix_table_F1.csv: ", nrow(f1), " rows")
}

message("Done.")
