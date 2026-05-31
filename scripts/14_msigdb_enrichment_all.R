#!/usr/bin/env Rscript
# 14_msigdb_enrichment_all.R
# MSigDB gene set enrichment for all 5 contrasts + ASO-specific DMRs
# Tests neural, synaptic, chromatin, and splicing gene sets
# Two independent enrichment methods: clusterProfiler (07_dmr_annotate.R)
# and MSigDB (this script) converging on same pathways = stronger evidence

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

setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/dmr_annotation/msigdb_v2"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

GROUP_COLS <- c(
  ASO_CTRL     = "#2E9B6F",
  ASO_VPA      = "#D94F3D",
  Scramble_VPA = "#F0A500",
  Scramble_CTRL= "#1D6FA4",
  ASO_specific = "#8E44AD"
)

# ── CONTRASTS TO ANALYSE ─────────────────────────────────────
CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       rds="results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds",
       label="ASO effect (nusinersen vs scramble CTRL)"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       rds="results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds",
       label="VPA effect (HDAC inhibitor vs CTRL)"),
  list(name="ASO_VPA_vs_Scramble_CTRL",
       rds="results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds",
       label="Combination vs baseline"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       rds="results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds",
       label="VPA effect on ASO background"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       rds="results/dmr/dmr_ASO_VPA_vs_Scramble_VPA.rds",
       label="ASO effect on VPA background"),
  list(name="ASO_specific",
       rds="results/dmr/dmr_ASO_specific.rds",
       label="ASO-specific DMRs (151 high-confidence)")
)

# ── GENE ANNOTATION FUNCTION ─────────────────────────────────
annotate_dmr <- function(rds_file, direction=NULL) {
  gr <- readRDS(rds_file)
  gr <- gr[mcols(gr)$cytosinesCount >= 6]
  if (!is.null(direction)) {
    if (direction == "hypo")
      gr <- gr[mcols(gr)$regionType == "gain"]
    else if (direction == "hyper")
      gr <- gr[mcols(gr)$regionType == "loss"]
  }
  if (length(gr) == 0) return(character(0))
  message("  Annotating ", length(gr), " DMRs...")
  anno <- annotatePeak(gr,
    TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
    tssRegion=c(-3000, 3000), verbose=FALSE)
  genes <- unique(as.character(anno@anno$geneId))
  genes[!is.na(genes) & genes != ""]
}

# ── MSIGDB GENE SETS ─────────────────────────────────────────
message("Loading MSigDB gene sets...")

# Neural / synaptic gene sets (C5 GO Biological Process)
neural_sets <- msigdbr(species="Homo sapiens", category="C5",
                       subcategory="GO:BP") |>
  dplyr::filter(grepl(
    "NEURON|SYNAP|AXON|DENDRIT|NEURAL|BRAIN|RETINAL|GANGLION|MOTOR_NEURON",
    gs_name, ignore.case=TRUE))

# Chromatin / epigenetic gene sets
chromatin_sets <- msigdbr(species="Homo sapiens", category="C5",
                          subcategory="GO:BP") |>
  dplyr::filter(grepl(
    "CHROMATIN|HISTONE|METHYLAT|ACETYL|EPIGENET|HETEROCHROMATIN",
    gs_name, ignore.case=TRUE))

# Splicing gene sets
splicing_sets <- msigdbr(species="Homo sapiens", category="C5",
                         subcategory="GO:BP") |>
  dplyr::filter(grepl(
    "SPLICING|SPLICEOSOM|MRNA_PROCESS|RNA_SPLICE",
    gs_name, ignore.case=TRUE))

# Reactome gene sets (C2)
reactome_sets <- msigdbr(species="Homo sapiens", category="C2",
                         subcategory="CP:REACTOME") |>
  dplyr::filter(grepl(
    "NEURON|SYNAP|AXON|NEURAL|NERVOUS",
    gs_name, ignore.case=TRUE))

make_term2gene <- function(sets) {
  data.frame(term=sets$gs_name, gene=sets$entrez_gene)
}

# ── RUN ENRICHMENT ───────────────────────────────────────────
run_msigdb <- function(genes, term2gene, min_genes=3) {
  if (length(genes) < min_genes) return(NULL)
  tryCatch(
    enricher(genes, TERM2GENE=term2gene,
             pAdjustMethod="BH", pvalueCutoff=0.05,
             minGSSize=5, maxGSSize=500),
    error=function(e) { message("  Error: ", e$message); NULL }
  )
}

save_result <- function(res, tag, contrast_name) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message("  No significant hits for: ", tag)
    return(invisible(NULL))
  }
  df <- as.data.frame(res)
  f <- file.path(OUT_DIR, paste0(contrast_name, "_", tag, ".tsv"))
  write.table(df, f, sep="\t", quote=FALSE, row.names=FALSE)
  message("  Saved: ", basename(f), " (", nrow(df), " terms)")
  df
}

make_dotplot <- function(res, title, col="#1B4F8A", top_n=15) {
  if (is.null(res)) return(NULL)
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  df <- head(df[order(df$p.adjust), ], top_n)
  df$Description <- gsub("GOBP_|REACTOME_", "", df$Description)
  df$Description <- gsub("_", " ", df$Description)
  df$Description <- tolower(df$Description)
  df$GeneRatio_num <- sapply(df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  ggplot(df, aes(x=GeneRatio_num, y=reorder(Description, -log10(p.adjust)),
                 size=Count, colour=p.adjust)) +
    geom_point() +
    scale_colour_gradient(low=col, high="grey70", name="p.adj") +
    scale_size_continuous(name="Gene count", range=c(2,8)) +
    labs(title=title, x="Gene ratio", y=NULL) +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold", size=10),
          axis.text.y=element_text(size=8))
}

# ── MAIN LOOP ────────────────────────────────────────────────
all_results <- list()

for (ct in CONTRASTS) {
  message("\n=== ", ct$name, " ===")

  # For ASO-specific use all DMRs; for contrasts split hyper/hypo
  if (ct$name == "ASO_specific") {
    genes_all <- annotate_dmr(ct$rds)
    directions <- list(list(dir=NULL, tag="all", genes=genes_all))
  } else {
    genes_hypo  <- annotate_dmr(ct$rds, "hypo")
    genes_hyper <- annotate_dmr(ct$rds, "hyper")
    directions  <- list(
      list(dir="hypo",  tag="hypo",  genes=genes_hypo),
      list(dir="hyper", tag="hyper", genes=genes_hyper)
    )
  }

  plots <- list()
  for (d in directions) {
    message(" Direction: ", d$tag, " (", length(d$genes), " genes)")
    if (length(d$genes) < 3) next

    res_neural    <- run_msigdb(d$genes, make_term2gene(neural_sets))
    res_chromatin <- run_msigdb(d$genes, make_term2gene(chromatin_sets))
    res_splicing  <- run_msigdb(d$genes, make_term2gene(splicing_sets))
    res_reactome  <- run_msigdb(d$genes, make_term2gene(reactome_sets))

    save_result(res_neural,    paste0("neural_",    d$tag), ct$name)
    save_result(res_chromatin, paste0("chromatin_", d$tag), ct$name)
    save_result(res_splicing,  paste0("splicing_",  d$tag), ct$name)
    save_result(res_reactome,  paste0("reactome_",  d$tag), ct$name)

    p_neural    <- make_dotplot(res_neural,
      paste0(ct$name, " — neural (", d$tag, ")"))
    p_chromatin <- make_dotplot(res_chromatin,
      paste0(ct$name, " — chromatin (", d$tag, ")"), col="#8E44AD")
    p_splicing  <- make_dotplot(res_splicing,
      paste0(ct$name, " — splicing (", d$tag, ")"), col="#D94F3D")
    p_reactome  <- make_dotplot(res_reactome,
      paste0(ct$name, " — Reactome neural (", d$tag, ")"), col="#2E9B6F")

    # Combined plot
    valid_plots <- Filter(Negate(is.null),
                          list(p_neural, p_reactome, p_chromatin, p_splicing))
    if (length(valid_plots) > 0) {
      combined <- wrap_plots(valid_plots, ncol=2) +
        plot_annotation(
          title=paste0("MSigDB enrichment: ", ct$name, " (", d$tag, ")"),
          subtitle=ct$label,
          theme=theme(plot.title=element_text(face="bold", size=12))
        )
      out_pdf <- file.path(OUT_DIR,
        paste0(ct$name, "_msigdb_", d$tag, "_combined.pdf"))
      ggsave(out_pdf, combined, width=14, height=10)
      message("  Saved: ", basename(out_pdf))
    } else {
      # Save negative result placeholder
      p_null <- ggplot() +
        annotate("text", x=0.5, y=0.5, size=5, colour="grey50",
          label=paste0("No significant enrichment detected\n",
                       "in MSigDB gene sets\n",
                       "(all p.adj > 0.05)\n",
                       ct$name, " (", d$tag, ")")) +
        theme_void() +
        labs(title=paste0("MSigDB: ", ct$name, " (", d$tag, ")"))
      out_pdf <- file.path(OUT_DIR,
        paste0(ct$name, "_msigdb_", d$tag, "_combined.pdf"))
      ggsave(out_pdf, p_null, width=8, height=5)
      message("  Saved: ", basename(out_pdf), " (negative result)")
    }
  }
}

message("\nAll done. Results in: ", OUT_DIR)
list.files(OUT_DIR, pattern="\\.pdf$")
