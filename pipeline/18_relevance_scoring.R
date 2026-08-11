#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/pairwise_context_scan"


aso_df <- read.csv(file.path(OUT, "ASO_context_dependent_top200.csv"))
vpa_df <- read.csv(file.path(OUT, "VPA_context_dependent_top200.csv"))



# 1. Known SMA modifiers and motor neuron genes
sma_relevant <- c(
  # Direct SMA
  "SMN1","SMN2","NAIP","SERF1A",
  # SMA modifiers (published)
  "NCALD","PLS3","NRXN2","CPNE1","TLL2","RHOB",
  # Motor neuron disease genes
  "SOD1","TDP43","TARDBP","FUS","C9orf72","ALS2","SETX","VAPB",
  "DCTN1","CHMP2B","UBQLN2","SQSTM1","VCP","MATR3","HNRNPA1",
  # SMA pathway - splicing
  "HNRNPA2B1","HNRNPD","PTBP1","PTBP2","RBM45","SRSF","NOVA1","NOVA2",
  # Neuromuscular junction
  "RAPSN","DOK7","MUSK","AGRN","LRP4","CHRNA1","CHRNB1","CHRND","CHRNE",
  # Motor neuron survival/development
  "BDNF","GDNF","CNTF","IGF1","NT3","VEGF","HGF",
  "ISL1","MNX1","CHAT","SLC18A3","SLC5A7",
  # Axon guidance (key for motor neurons)
  "SEMA3A","SEMA3C","SEMA3F","NRPN1","NRPN2","PLXNA1","PLXNA2",
  "EPHB1","EPHB2","EPHA4","EFNB1","EFNB2","ROBO1","ROBO2","SLIT1",
  "UNC5A","UNC5B","UNC5C","DCC","NTN1",
  # Chromatin/epigenetic (VPA targets)
  "HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8",
  "KAT2A","KAT2B","EP300","CREBBP","BRD4","DNMT1","DNMT3A","DNMT3B",
  "EZH2","SUZ12","EED","KDM1A","KDM5C","KDM6A",
  # Calcium/kinase signalling (neural)
  "CAMK2A","CAMK2B","CAMK4","DYRK1A","CDK5","GSK3B",
  # RNA processing
  "SMN1","GEMIN2","GEMIN3","GEMIN4","GEMIN5","UNRIP","WDR77",
  # Ubiquitin/proteasome (USP7 etc)
  "USP7","USP27X","UBE3A","TRIM32","TRIM71","FBXO32","MURF1",
  # Transcription factors relevant to SMA/motor neurons
  "IRF8","PAX5","SOX11","ETV1","ETV4","ETV5","PEA3"
)

# 2. GO term based scoring — neural/motor neuron terms
neural_go_terms <- c(
  "neuromuscular junction", "motor neuron", "axon guidance",
  "neuron differentiation", "synaptic transmission", "axonogenesis",
  "neuron projection", "dendritic spine", "calcium signaling",
  "RNA splicing", "chromatin remodeling", "histone modification",
  "ubiquitin", "proteasome"
)


score_genes <- function(df, label) {
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  df <- df[!grepl("^LOC|^LINC|^MIR|^SNOR|-AS[0-9]|-DT$", df$SYMBOL), ]

  df$score <- 0

  # Score 1: known SMA/motor neuron relevant gene (+3)
  df$score <- df$score + ifelse(df$SYMBOL %in% sma_relevant, 3, 0)
  df$sma_relevant <- df$SYMBOL %in% sma_relevant

  # Score 2: promoter location (+2, more likely functional)
  df$score <- df$score + ifelse(grepl("Promoter", df$annotation), 2, 0)

  # Score 3: effect size rank (top quartile = +2)
  df$effect_rank <- rank(-df$meth_diff) / nrow(df)
  df$score <- df$score + ifelse(df$effect_rank <= 0.25, 2, 0)
  df$score <- df$score + ifelse(df$effect_rank <= 0.10, 1, 0)

  # Score 4: not intergenic (+1)
  df$score <- df$score + ifelse(!grepl("Intergenic", df$annotation), 1, 0)

  df$group <- label
  df[order(df$score, df$meth_diff, decreasing=TRUE), ]
}

message("Scoring ASO context-dependent genes...")
aso_scored <- score_genes(aso_df, "ASO_context_dependent")

message("Scoring VPA context-dependent genes...")
vpa_scored <- score_genes(vpa_df, "VPA_context_dependent")


write.csv(aso_scored, file.path(OUT, "ASO_context_dependent_scored.csv"),
          row.names=FALSE)
write.csv(vpa_scored, file.path(OUT, "VPA_context_dependent_scored.csv"),
          row.names=FALSE)


cat("\n=== TOP ASO CONTEXT-DEPENDENT GENES (by relevance score) ===\n")
top_aso <- aso_scored[aso_scored$score >= 2,
  c("SYMBOL","seqnames","meth_diff","annotation","score","sma_relevant")]
print(head(top_aso, 30))

cat("\n=== TOP VPA CONTEXT-DEPENDENT GENES (by relevance score) ===\n")
top_vpa <- vpa_scored[vpa_scored$score >= 2,
  c("SYMBOL","seqnames","meth_diff","annotation","score","sma_relevant")]
print(head(top_vpa, 30))


run_go <- function(df, label) {
  genes <- unique(df$SYMBOL[df$score >= 2])
  if (length(genes) < 5) {
    message("Too few high-scoring genes for GO: ", label)
    return(NULL)
  }
  ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  go  <- enrichGO(gene=ids$ENTREZID, OrgDb=org.Hs.eg.db,
                  ont="BP", pAdjustMethod="BH",
                  pvalueCutoff=0.05, readable=TRUE)
  if (is.null(go) || nrow(as.data.frame(go)) == 0) {
    message("No GO terms for: ", label); return(NULL)
  }
  message(label, ": ", nrow(as.data.frame(go)), " GO terms")
  go
}

message("\nRunning GO enrichment on high-scoring genes...")
go_aso <- run_go(aso_scored, "ASO_context")
go_vpa <- run_go(vpa_scored, "VPA_context")


make_bubble <- function(df, title, col) {
  top <- head(df[df$score >= 2, ], 25)
  top <- top[!duplicated(top$SYMBOL), ]
  if (nrow(top) == 0) return(NULL)
  top$SYMBOL <- factor(top$SYMBOL,
                       levels=top$SYMBOL[order(top$score, top$meth_diff)])
  ggplot(top, aes(x=meth_diff, y=SYMBOL, size=score,
                  colour=sma_relevant)) +
    geom_point(alpha=0.8) +
    scale_colour_manual(values=c("TRUE"="#E31A1C", "FALSE"=col),
                        labels=c("TRUE"="SMA/neural relevant",
                                 "FALSE"="Other"),
                        name=NULL) +
    scale_size_continuous(range=c(3,8), name="Relevance\nscore") +
    labs(title=title,
         x="|Methylation difference|", y=NULL) +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold", size=10),
          legend.position="right")
}

p1 <- make_bubble(aso_scored, "ASO context-dependent top hits", "#1F3A5F")
p2 <- make_bubble(vpa_scored, "VPA context-dependent top hits", "#C0392B")

if (!is.null(p1) && !is.null(p2)) {
  combined <- (p1 | p2) +
    plot_annotation(
      title="Context-dependent DMR candidates — relevance scored",
      caption="Red = known SMA/motor neuron/chromatin relevant gene. Size = relevance score.",
      theme=theme(plot.title=element_text(face="bold", size=13))
    )
  ggsave(file.path(OUT, "context_candidates_bubble.pdf"),
         combined, width=16, height=10, device=cairo_pdf)
  ggsave(file.path(OUT, "context_candidates_bubble.png"),
         combined, width=16, height=10, dpi=150)
  message("Saved bubble chart")
}


if (!is.null(go_aso)) {
  pdf(file.path(OUT, "GO_ASO_context_highscore.pdf"), width=10, height=8)
  print(dotplot(go_aso, showCategory=15,
                title="GO BP: high-scoring ASO context-dependent genes"))
  dev.off()
}
if (!is.null(go_vpa)) {
  pdf(file.path(OUT, "GO_VPA_context_highscore.pdf"), width=10, height=8)
  print(dotplot(go_vpa, showCategory=15,
                title="GO BP: high-scoring VPA context-dependent genes"))
  dev.off()
}

message("\nAll done. Outputs in: ", OUT)
