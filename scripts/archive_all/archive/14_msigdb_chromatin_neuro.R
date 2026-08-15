suppressPackageStartupMessages({
  library(GenomicRanges)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(patchwork)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})
.libPaths(c("~/R/library", .libPaths()))
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/dmr_annotation/msigdb"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

annotate_dmr <- function(rds_file, filter_fn=NULL) {
  gr <- readRDS(rds_file)
  gr <- gr[mcols(gr)$cytosinesCount >= 6]
  if (!is.null(filter_fn)) gr <- filter_fn(gr)
  message("  Annotating ", length(gr), " DMRs...")
  anno <- annotatePeak(gr,
    TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
    tssRegion=c(-3000, 3000), verbose=FALSE)
  genes <- unique(as.character(anno@anno$geneId))
  genes[!is.na(genes) & genes != ""]
}

message("Loading and annotating DMR gene lists...")
# VPA hypomethylated DMRs — expect HDAC/chromatin enrichment
genes_vpa_hypo <- annotate_dmr(
  "results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds",
  function(gr) gr[mcols(gr)$direction == 1])

# ASO hypermethylated DMRs — expect neurodegeneration/motor neuron enrichment
genes_aso_hyper <- annotate_dmr(
  "results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds",
  function(gr) gr[mcols(gr)$direction == 1])

# ASO-specific DMRs — all directions
genes_specific <- annotate_dmr("results/dmr/dmr_ASO_specific.rds")

message("  VPA hypo genes: ",    length(genes_vpa_hypo))
message("  ASO hyper genes: ",   length(genes_aso_hyper))
message("  ASO-specific genes: ", length(genes_specific))

message("Fetching MSigDB gene sets...")

# Analysis 1: VPA vs HDAC/chromatin remodelling
# C5 GO:BP terms related to histone modification and chromatin
msig_c5 <- msigdbr(species="Homo sapiens",
                   collection="C5", subcollection="GO:BP")
msig_h  <- msigdbr(species="Homo sapiens", collection="H")

chromatin_terms <- c(
  "HISTONE","CHROMATIN","HDAC","DEACETYL","ACETYL",
  "METHYLAT","EPIGENET","NUCLEOSOM","HETEROCHROMATIN",
  "REMODEL","DNMT","POLYCOMB","TRITHORAX"
)
chromatin_sets <- rbind(
  msig_c5[grepl(paste(chromatin_terms, collapse="|"),
                toupper(msig_c5$gs_name)), ],
  msig_h[grepl(paste(chromatin_terms, collapse="|"),
               toupper(msig_h$gs_name)), ]
)
message("  Chromatin gene sets: ", length(unique(chromatin_sets$gs_name)))

# Analysis 2: ASO hyper DMRs vs neurodegeneration/motor neuron
msig_c2 <- msigdbr(species="Homo sapiens",
                   collection="C2", subcollection="CP:REACTOME")
msig_c2_cgp <- msigdbr(species="Homo sapiens",
                        collection="C2", subcollection="CGP")

neuro_terms <- c(
  "NEURON","MOTOR_NEURON","AXON","SYNAP","NEURODEGENERAT",
  "AMYOTROPH","SPINAL","SMN","SMA","MUSCULAR_ATROPHY",
  "DENDRIT","NEURODEVELOP","AXONOGENESIS"
)
neuro_sets <- rbind(
  msig_c2[grepl(paste(neuro_terms, collapse="|"),
               toupper(msig_c2$gs_name)), ],
  msig_c2_cgp[grepl(paste(neuro_terms, collapse="|"),
                    toupper(msig_c2_cgp$gs_name)), ],
  msig_c5[grepl(paste(neuro_terms, collapse="|"),
               toupper(msig_c5$gs_name)), ]
)
message("  Neuro gene sets: ", length(unique(neuro_sets$gs_name)))

make_term2gene <- function(sets) {
  t2g <- sets[, c("gs_name", "ncbi_gene")]
  t2g$ncbi_gene <- as.character(t2g$ncbi_gene)
  t2g
}

run_enrich <- function(genes, term2gene, label) {
  message("  Running enrichment: ", label, " (n=", length(genes), " genes)")
  res <- enricher(
    gene          = genes,
    TERM2GENE     = term2gene,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    minGSSize     = 5,
    maxGSSize     = 500
  )
  n_sig <- if (!is.null(res)) sum(as.data.frame(res)$p.adjust <= 0.05) else 0
  message("    Significant terms: ", n_sig)
  res
}

t2g_chromatin <- make_term2gene(chromatin_sets)
t2g_neuro     <- make_term2gene(neuro_sets)

res_vpa_chromatin  <- run_enrich(genes_vpa_hypo,  t2g_chromatin, "VPA hypo vs chromatin")
res_aso_neuro      <- run_enrich(genes_aso_hyper,  t2g_neuro,     "ASO hyper vs neuro")
res_specific_neuro <- run_enrich(genes_specific,   t2g_neuro,     "ASO-specific vs neuro")

# Save TSVs
save_tsv <- function(res, prefix) {
  df <- if (!is.null(res)) as.data.frame(res) else data.frame()
  write.table(df,
    file.path(OUT_DIR, paste0(prefix, "_msigdb.tsv")),
    sep="\t", quote=FALSE, row.names=FALSE)
  df
}

df_vpa  <- save_tsv(res_vpa_chromatin,  "VPA_hypo_chromatin")
df_aso  <- save_tsv(res_aso_neuro,      "ASO_hyper_neuro")
df_spec <- save_tsv(res_specific_neuro, "ASO_specific_neuro")

# Plot function — always produces a visible figure
make_plot <- function(res, label, subtitle) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    return(
      ggplot(data.frame(x=0.5, y=0.5), aes(x, y)) +
        annotate("label", x=0.5, y=0.5, size=4.5, fill="grey95",
                 color="grey30", label.size=0.4,
                 label="No significant enrichment detected") +
        theme_classic() +
        theme(axis.text=element_blank(), axis.ticks=element_blank(),
              axis.line=element_blank(), axis.title=element_blank()) +
        labs(title=label, subtitle=subtitle)
    )
  }
  dotplot(res, showCategory=20, font.size=9) +
    labs(title=label, subtitle=subtitle) +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold", size=11),
          plot.subtitle=element_text(size=9, color="grey40"))
}

p1 <- make_plot(res_vpa_chromatin,
                "VPA hypomethylated DMRs — chromatin/HDAC gene sets",
                "Scramble_VPA vs Scramble_CTRL · C5 GO:BP + Hallmark")

p2 <- make_plot(res_aso_neuro,
                "ASO hypermethylated DMRs — neurodegeneration/motor neuron gene sets",
                "ASO_CTRL vs Scramble_CTRL · C2 REACTOME + CGP + C5 GO:BP")

p3 <- make_plot(res_specific_neuro,
                "ASO-specific DMRs — neurodegeneration/motor neuron gene sets",
                "151 ASO-specific DMRs · C2 REACTOME + CGP + C5 GO:BP")

ggsave(file.path(OUT_DIR, "msigdb_VPA_chromatin.pdf"),
       p1, width=11, height=8)
ggsave(file.path(OUT_DIR, "msigdb_ASO_neuro.pdf"),
       p2, width=11, height=8)
ggsave(file.path(OUT_DIR, "msigdb_ASO_specific_neuro.pdf"),
       p3, width=11, height=8)

message("\nSummary:")
for (nm in c("VPA chromatin", "ASO neuro", "ASO-specific neuro")) {
  res <- list(res_vpa_chromatin, res_aso_neuro, res_specific_neuro)[[
    which(c("VPA chromatin","ASO neuro","ASO-specific neuro") == nm)]]
  n <- if (!is.null(res)) sum(as.data.frame(res)$p.adjust <= 0.05) else 0
  message("  ", nm, ": ", n, " significant terms")
}
message("Done. Outputs in: ", OUT_DIR)
