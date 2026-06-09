.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/dmr_annotation"

# ── ITEM 10: Chr13 hotspot gene annotation ────────────────────────────────────
message('\n=== ITEM 10: Chromosome 13 DMR hotspot annotation ===')

# Chr13 hotspot region from meeting: ~60-80 Mb
CHR13_START <- 60000000
CHR13_END   <- 80000000

# Load ASO_VPA vs Scramble_CTRL DMRs (largest contrast)
dmr_file <- 'results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds'
if (file.exists(dmr_file)) {
  dmrs <- readRDS(dmr_file)
  message('Loaded DMRs: ', length(dmrs))

  # Filter to chr13 hotspot
  chr13_dmrs <- dmrs[seqnames(dmrs)=='chr13' &
                     start(dmrs) >= CHR13_START &
                     end(dmrs)   <= CHR13_END]
  message('Chr13 hotspot DMRs (60-80Mb): ', length(chr13_dmrs))

  if (length(chr13_dmrs) > 0) {
    # Annotate
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    anno <- annotatePeak(chr13_dmrs, tssRegion=c(-2000,2000),
      TxDb=txdb, annoDb='org.Hs.eg.db', verbose=FALSE)
    anno_df <- as.data.frame(anno)

    # Get unique genes
    genes <- unique(anno_df$SYMBOL[!is.na(anno_df$SYMBOL)])
    message('Unique genes in chr13 hotspot: ', length(genes))
    message('Top genes: ', paste(head(genes, 20), collapse=', '))

    write.csv(anno_df,
      'results/dmr_annotation/chr13_hotspot_annotated.csv',
      row.names=FALSE)

    # GO enrichment
    gene_ids <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
    if (length(gene_ids) >= 10) {
      go_res <- enrichGO(gene=gene_ids, OrgDb=org.Hs.eg.db,
        ont='BP', pAdjustMethod='BH',
        pvalueCutoff=0.05, readable=TRUE)

      if (!is.null(go_res) && nrow(go_res@result[go_res@result$p.adjust<0.05,]) > 0) {
        message('\nSignificant GO terms in chr13 hotspot:')
        print(head(go_res@result[go_res@result$p.adjust<0.05,
          c('Description','p.adjust','Count')], 10))
        write.csv(go_res@result,
          'results/dmr_annotation/chr13_hotspot_GO.csv',
          row.names=FALSE)

        pdf('results/dmr_annotation/chr13_hotspot_GO_dotplot.pdf',
            width=10, height=8)
        print(dotplot(go_res, showCategory=15,
          title='GO BP — Chr13 DMR hotspot (60-80Mb)\nASO_VPA vs Scramble_CTRL'))
        dev.off()
      } else {
        message('No significant GO terms in chr13 hotspot')
      }
    }

    # Also check chromatin/heterochromatin genes specifically
    message('\n--- Checking for heterochromatin/repeat genes ---')
    hetero_genes <- anno_df$SYMBOL[grep('KRAB|ZNF|LMNB|SATB|HP1|CBX|HDAC|DNMT|H3K|KDM|KMT',
      anno_df$SYMBOL, ignore.case=TRUE)]
    if (length(hetero_genes) > 0) {
      message('Chromatin-related genes: ', paste(unique(hetero_genes), collapse=', '))
    }

    # Region type summary
    message('\nAnnotation summary:')
    print(table(anno_df$annotation))
  }
} else {
  message('DMR file not found: ', dmr_file)
  # Try alternative
  dmr_files <- list.files('results/dmr', pattern='*.rds', full.names=TRUE)
  message('Available DMR files: ', paste(dmr_files, collapse=', '))
}

message('\nDone.')


dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

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
       label="ASO-specific DMRs")
)

# annotate DMRs and return Entrez gene IDs
annotate_dmr <- function(rds_file, direction=NULL) {
  gr <- readRDS(rds_file)
  gr <- gr[mcols(gr)$cytosinesCount >= 6]
  if (!is.null(direction)) {
    gr <- gr[mcols(gr)$regionType == if(direction=="hypo") "gain" else "loss"]
  }
  if (length(gr) == 0) return(character(0))
  message("  annotating ", length(gr), " DMRs...")
  anno  <- annotatePeak(gr, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                        tssRegion=c(-3000,3000), verbose=FALSE)
  genes <- unique(as.character(anno@anno$geneId))
  genes[!is.na(genes) & genes != ""]
}

# load MSigDB gene sets
message("loading MSigDB gene sets...")
neural_sets    <- msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP") |>
  dplyr::filter(grepl("NEURON|SYNAP|AXON|DENDRIT|NEURAL|BRAIN|MOTOR_NEURON",
                      gs_name, ignore.case=TRUE))
chromatin_sets <- msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP") |>
  dplyr::filter(grepl("CHROMATIN|HISTONE|METHYLAT|ACETYL|EPIGENET|HETEROCHROMATIN",
                      gs_name, ignore.case=TRUE))
splicing_sets  <- msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP") |>
  dplyr::filter(grepl("SPLICING|SPLICEOSOM|MRNA_PROCESS|RNA_SPLICE",
                      gs_name, ignore.case=TRUE))
reactome_sets  <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME") |>
  dplyr::filter(grepl("NEURON|SYNAP|AXON|NEURAL|NERVOUS",
                      gs_name, ignore.case=TRUE))

make_term2gene <- function(sets) data.frame(term=sets$gs_name, gene=sets$entrez_gene)

run_msigdb <- function(genes, term2gene, min_genes=3) {
  if (length(genes) < min_genes) return(NULL)
  tryCatch(
    enricher(genes, TERM2GENE=term2gene, pAdjustMethod="BH",
             pvalueCutoff=0.05, minGSSize=5, maxGSSize=500),
    error=function(e) { message("  error: ", e$message); NULL }
  )
}

save_result <- function(res, tag, contrast_name) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message("  no hits: ", tag); return(invisible(NULL))
  }
  df <- as.data.frame(res)
  f  <- file.path(OUT_DIR, paste0(contrast_name, "_", tag, ".tsv"))
  write.table(df, f, sep="\t", quote=FALSE, row.names=FALSE)
  message("  saved: ", basename(f), " (", nrow(df), " terms)")
  df
}

make_dotplot <- function(res, title, col="#1B4F8A", top_n=15) {
  if (is.null(res)) return(NULL)
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  df <- head(df[order(df$p.adjust), ], top_n)
  df$Description    <- tolower(gsub("_", " ", gsub("GOBP_|REACTOME_", "", df$Description)))
  df$GeneRatio_num  <- sapply(df$GeneRatio, function(x) {
    p <- strsplit(x, "/")[[1]]; as.numeric(p[1])/as.numeric(p[2]) })
  ggplot(df, aes(x=GeneRatio_num, y=reorder(Description, -log10(p.adjust)),
                 size=Count, colour=p.adjust)) +
    geom_point() +
    scale_colour_gradient(low=col, high="grey70", name="p.adj") +
    scale_size_continuous(name="genes", range=c(2,8)) +
    labs(title=title, x="gene ratio", y=NULL) +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold", size=10),
          axis.text.y=element_text(size=8))
}

# main loop
for (ct in CONTRASTS) {
  message("\n", ct$name)
  directions <- if (ct$name == "ASO_specific") {
    list(list(tag="all", genes=annotate_dmr(ct$rds)))
  } else {
    list(list(tag="hypo",  genes=annotate_dmr(ct$rds, "hypo")),
         list(tag="hyper", genes=annotate_dmr(ct$rds, "hyper")))
  }

  for (d in directions) {
    message("  direction: ", d$tag, " (", length(d$genes), " genes)")
    if (length(d$genes) < 3) next

    res_neural    <- run_msigdb(d$genes, make_term2gene(neural_sets))
    res_chromatin <- run_msigdb(d$genes, make_term2gene(chromatin_sets))
    res_splicing  <- run_msigdb(d$genes, make_term2gene(splicing_sets))
    res_reactome  <- run_msigdb(d$genes, make_term2gene(reactome_sets))

    save_result(res_neural,    paste0("neural_",    d$tag), ct$name)
    save_result(res_chromatin, paste0("chromatin_", d$tag), ct$name)
    save_result(res_splicing,  paste0("splicing_",  d$tag), ct$name)
    save_result(res_reactome,  paste0("reactome_",  d$tag), ct$name)

    valid_plots <- Filter(Negate(is.null), list(
      make_dotplot(res_neural,    paste0("neural (", d$tag, ")")),
      make_dotplot(res_reactome,  paste0("reactome (", d$tag, ")"), col="#2E9B6F"),
      make_dotplot(res_chromatin, paste0("chromatin (", d$tag, ")"), col="#8E44AD"),
      make_dotplot(res_splicing,  paste0("splicing (", d$tag, ")"),  col="#D94F3D")
    ))

    out_pdf <- file.path(OUT_DIR, paste0(ct$name, "_msigdb_", d$tag, "_combined.pdf"))

    if (length(valid_plots) > 0) {
      combined <- wrap_plots(valid_plots, ncol=2) +
        plot_annotation(title=paste0("MSigDB: ", ct$name, " (", d$tag, ")"),
                        subtitle=ct$label,
                        theme=theme(plot.title=element_text(face="bold", size=12)))
      ggsave(out_pdf, combined, width=14, height=10)
    } else {
      p_null <- ggplot() +
        annotate("text", x=0.5, y=0.5, size=5, colour="grey50",
                 label=paste0("no significant enrichment\n", ct$name, " (", d$tag, ")")) +
        theme_void()
      ggsave(out_pdf, p_null, width=8, height=5)
    }
    message("  saved: ", basename(out_pdf))
  }
}
message("\ndone. results in: ", OUT_DIR)
