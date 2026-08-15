#!/usr/bin/env Rscript
# 17_pairwise_candidates_unified.R
# Pairwise candidate gene selection for SMA epigenomics project
#
# Four pairwise contrasts (locked params: binSize=300, minDiff=0.20, p<0.01):
#   ASO alone:  ASO_CTRL vs Scramble_CTRL
#   VPA alone:  Scramble_VPA vs Scramble_CTRL
#   ASO in VPA: ASO_VPA vs Scramble_VPA
#   VPA in ASO: ASO_VPA vs ASO_CTRL
#
# Candidate categories:
#   synergy     -- gene has DMR in BOTH combination contrasts,
#                  absent from its corresponding single-drug contrast
#   ASO_background_specific -- DMR in ASO_in_VPA only (not VPA_in_ASO, not ASO_alone)
#   VPA_background_specific -- DMR in VPA_in_ASO only (not ASO_in_VPA, not VPA_alone)
#
# Input:  results/dmr/dmr_<contrast>.rds (one per contrast)
# Output: results/pairwise_context_scan/
#
# Usage: Rscript 17_pairwise_candidates_unified.R
#        or submit via sbatch (64G, 2h)

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(dplyr)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/pairwise_context_scan"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 1. Load DMR sets
message("Loading DMR RDS files...")
aso_alone  <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
vpa_alone  <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
aso_in_vpa <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_VPA.rds")
vpa_in_aso <- readRDS("results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds")

cat("DMR counts:\n")
cat("  ASO alone:", length(aso_alone), "\n")
cat("  VPA alone:", length(vpa_alone), "\n")
cat("  ASO in VPA:", length(aso_in_vpa), "\n")
cat("  VPA in ASO:", length(vpa_in_aso), "\n")

# 2. Annotate with ChIPseeker
#    tssRegion = +/-2kb around TSS
#    proportion1 = treatment, proportion2 = reference (see 02_dmr_calling.R)
#    meth_diff > 0 means treatment hypermethylated (regionType "loss")
annotate_dmrs <- function(gr, label) {
  anno <- annotatePeak(gr, tssRegion=c(-2000, 2000),
                       TxDb=txdb, annoDb="org.Hs.eg.db", verbose=FALSE)
  df <- as.data.frame(anno)
  df$contrast  <- label
  df$meth_diff <- df$proportion1 - df$proportion2
  df
}

message("Annotating...")
ann_aso_alone  <- annotate_dmrs(aso_alone,  "ASO_alone")
ann_vpa_alone  <- annotate_dmrs(vpa_alone,  "VPA_alone")
ann_aso_in_vpa <- annotate_dmrs(aso_in_vpa, "ASO_in_VPA")
ann_vpa_in_aso <- annotate_dmrs(vpa_in_aso, "VPA_in_ASO")

write.csv(ann_aso_alone,  file.path(OUT, "annotated_ASO_alone.csv"),  row.names=FALSE)
write.csv(ann_vpa_alone,  file.path(OUT, "annotated_VPA_alone.csv"),  row.names=FALSE)
write.csv(ann_aso_in_vpa, file.path(OUT, "annotated_ASO_in_VPA.csv"), row.names=FALSE)
write.csv(ann_vpa_in_aso, file.path(OUT, "annotated_VPA_in_ASO.csv"), row.names=FALSE)

# 3. Gene-level sets
#    Filter out uncharacterised loci (LOC*, LINC*, MIR*, antisense)
clean_genes <- function(df) {
  unique(df$SYMBOL[!is.na(df$SYMBOL) & df$SYMBOL != "" &
                   !grepl("^LOC|^LINC|^MIR|^SNOR|-AS[0-9]$|-DT$", df$SYMBOL)])
}

genes_aso_alone  <- clean_genes(ann_aso_alone)
genes_vpa_alone  <- clean_genes(ann_vpa_alone)
genes_aso_in_vpa <- clean_genes(ann_aso_in_vpa)
genes_vpa_in_aso <- clean_genes(ann_vpa_in_aso)

# 4. Candidate categories
#    context-dependent: present in combination contrast, absent in single-drug
#    synergy:           context-dependent in both combination contrasts
#    restricted:        context-dependent in only one combination contrast
aso_context_dep <- setdiff(genes_aso_in_vpa, genes_aso_alone)
vpa_context_dep <- setdiff(genes_vpa_in_aso, genes_vpa_alone)
synergy         <- intersect(aso_context_dep, vpa_context_dep)
aso_restricted  <- setdiff(aso_context_dep, vpa_context_dep)
vpa_restricted  <- setdiff(vpa_context_dep, aso_context_dep)

cat("\nCandidate counts:\n")
cat("  ASO context-dependent:", length(aso_context_dep), "\n")
cat("  VPA context-dependent:", length(vpa_context_dep), "\n")
cat("  Pairwise synergy:", length(synergy), "\n")
cat("  ASO restricted:", length(aso_restricted), "\n")
cat("  VPA restricted:", length(vpa_restricted), "\n")

# 5. Relevance scoring (0-8 per DMR; take max across DMRs per gene)
#    +3 SMA/motor neuron relevant gene
#    +3 promoter; +2 UTR/exon; +1 intron/intergenic
#    +2 |meth_diff| >= 0.40; +1 >= 0.20
sma_relevant <- c(
  "SMN1","SMN2","NAIP","NCALD","PLS3",
  "SOD1","TARDBP","FUS","ALS2",
  "HNRNPA1","HNRNPA2B1","PTBP1","PTBP2","SRSF1","NOVA1","NOVA2",
  "CHAT","SLC18A3","ISL1","MNX1",
  "SEMA3A","EPHB1","EPHB2","EPHA4","ROBO1","ROBO2","NTN1","DCC",
  "BDNF","GDNF","CNTF","IGF1",
  "HDAC1","HDAC2","HDAC3","HDAC4","HDAC6","KDM1A","EZH2",
  "CAMK2A","CAMK2B","CAMK4","CDK5","GSK3B",
  "USP7","USP27X","UBE3A","TRIM32",
  "IRF8","PAX5","SOX11","ETV1","ETV4",
  "SHANK2","SHANK3","DLGAP1","DLG4","NRXN1","CNTNAP2"
)

score_gene <- function(sym, ann_df) {
  rows <- ann_df[ann_df$SYMBOL == sym & !is.na(ann_df$SYMBOL), ]
  if (nrow(rows) == 0) return(0)
  s <- 0
  s <- s + ifelse(any(grepl("Promoter", rows$annotation)), 3,
           ifelse(any(grepl("5' UTR|Exon", rows$annotation)), 2, 1))
  s <- s + ifelse(sym %in% sma_relevant, 3, 0)
  s <- s + ifelse(any(abs(rows$meth_diff) >= 0.40), 2,
           ifelse(any(abs(rows$meth_diff) >= 0.20), 1, 0))
  s
}

score_set <- function(genes, ann_df, label) {
  scores <- sapply(genes, function(g) score_gene(g, ann_df))
  df <- data.frame(SYMBOL=genes, score=scores, category=label,
                   stringsAsFactors=FALSE)
  df[order(df$score, decreasing=TRUE), ]
}

# synergy: score against both combination contrasts, take maximum
synergy_scores <- sapply(synergy, function(g) {
  max(score_gene(g, ann_aso_in_vpa), score_gene(g, ann_vpa_in_aso))
})
synergy_df     <- data.frame(SYMBOL=synergy, score=synergy_scores,
                              category="synergy", stringsAsFactors=FALSE)
synergy_df     <- synergy_df[order(synergy_df$score, decreasing=TRUE), ]
aso_rest_df    <- score_set(aso_restricted, ann_aso_in_vpa, "ASO_background_specific")
vpa_rest_df    <- score_set(vpa_restricted, ann_vpa_in_aso, "VPA_background_specific")

cat("\nTop synergy candidates (score >= 3):\n")
print(head(synergy_df[synergy_df$score >= 3, ], 20))
cat("\nTop ASO-restricted (score >= 3):\n")
print(head(aso_rest_df[aso_rest_df$score >= 3, ], 10))
cat("\nTop VPA-restricted (score >= 3):\n")
print(head(vpa_rest_df[vpa_rest_df$score >= 3, ], 10))

write.csv(synergy_df,    file.path(OUT, "synergy_candidates_pairwise.csv"),  row.names=FALSE)
write.csv(aso_rest_df,   file.path(OUT, "ASO_background_specific_pairwise.csv"),      row.names=FALSE)
write.csv(vpa_rest_df,   file.path(OUT, "VPA_background_specific_pairwise.csv"),      row.names=FALSE)

message("Done. Outputs: ", OUT)
