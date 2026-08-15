.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
source("scripts/pipeline/00_sma_palette.R")

# Filter to bona fide TFs using GO:0003700 (DNA-binding TF activity)
library(org.Hs.eg.db)
tf_entrez <- AnnotationDbi::select(org.Hs.eg.db,
  keys="GO:0003700", keytype="GOALL", columns="SYMBOL")$SYMBOL
tf_whitelist <- unique(toupper(tf_entrez[!is.na(tf_entrez)]))
message("TF whitelist: ", length(tf_whitelist), " genes")
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/tf_motif"
TOPN    <- 30

extract_gene_symbol <- function(target) {
  toks <- unlist(strsplit(target, "[_./-]"))
  toks <- toks[grepl("^[A-Za-z][A-Za-z0-9]*$", toks)]
  if (length(toks) == 0) return(NA_character_)
  toupper(toks[length(toks)])
}

CONTRASTS <- list(
  list(name="ASO",     csv=file.path(OUT_DIR,"pwmenrich_top_motifs_ASO.csv")),
  list(name="VPA",     csv=file.path(OUT_DIR,"pwmenrich_top_motifs_VPA.csv")),
  list(name="ASO_VPA", csv=file.path(OUT_DIR,"pwmenrich_top_motifs_ASO_VPA.csv"))
)

message("Initialising STRINGdb...")
string_db <- STRINGdb$new(version="12.0", species=9606,
                          score_threshold=400, input_directory="")

build_net <- function(name, csv_path) {
  message("\n=== ", name, " ===")
  if (!file.exists(csv_path)) { message("  missing: ", csv_path); return(NULL) }
  df <- read.csv(csv_path, stringsAsFactors=FALSE)
  df <- df[order(df$p.value),]
  df$gene <- sapply(df$target, extract_gene_symbol)
  df <- df[!is.na(df$gene) & nchar(df$gene)>=2,]
  df <- df[!duplicated(df$gene),]
  df <- df[df$gene %in% tf_whitelist,]
  message("  after TF filter: ", nrow(df), " genes")
  topdf <- head(df, TOPN)
  message("  genes: ", paste(topdf$gene, collapse=", "))
  mapped <- string_db$map(topdf, "gene", removeUnmappedRows=TRUE)
  message("  mapped: ", nrow(mapped), "/", nrow(topdf))
  if (nrow(mapped) < 3) return(NULL)
  ints <- string_db$get_interactions(mapped$STRING_id)
  if (nrow(ints) == 0) return(NULL)
  g <- graph_from_data_frame(ints[,c("from","to")], directed=FALSE,
       vertices=data.frame(name=mapped$STRING_id, gene=mapped$gene))
  g <- igraph::simplify(g)
  g <- delete_vertices(g, which(igraph::degree(g)==0))
  if (vcount(g) < 3) return(NULL)
  V(g)$degree <- igraph::degree(g)
  V(g)$label  <- V(g)$gene
  deg_df <- data.frame(gene=V(g)$gene, degree=V(g)$degree)
  deg_df <- deg_df[order(-deg_df$degree),]
  write.csv(deg_df, file.path(OUT_DIR, paste0("tf_hub_degree_",name,".csv")), row.names=FALSE)
  message("  hubs: ", paste(head(deg_df$gene,5), collapse=", "))
  list(graph=g, contrast=name)
}

render_panel <- function(net) {
  g <- net$graph
  set.seed(42)
  ggraph(g, layout="fr") +
    geom_edge_link(colour="grey70", alpha=0.45, edge_width=0.35) +
    geom_node_point(aes(size=degree, fill=degree), shape=21, colour="grey25", stroke=0.4) +
    geom_node_text(aes(label=label), size=3, colour="black", fontface="bold", repel=TRUE) +
    scale_fill_gradient(low="#FFE066", high="#C0392B", name="Degree") +
    scale_size_continuous(range=c(4,12), guide="none") +
    theme_void(base_size=12) +
    theme(plot.title=element_text(face="bold", hjust=0.5), legend.position="right") +
    labs(title=paste0("TF PPI hub network: ", net$contrast),
         subtitle="STRINGdb v12, top 30 PWMEnrich motifs, node=degree")
}

nets <- list()
for (ct in CONTRASTS) {
  net <- build_net(ct$name, ct$csv)
  if (!is.null(net)) {
    nets[[ct$name]] <- net
    p <- render_panel(net)
    ggsave(file.path(OUT_DIR, paste0("tf_hub_network_",ct$name,".pdf")),
           p, width=10, height=9, device=cairo_pdf)
    message("saved: tf_hub_network_", ct$name, ".pdf")
  }
}

if (length(nets) >= 2) {
  panels <- lapply(names(nets), function(nm) render_panel(nets[[nm]]))
  combined <- wrap_plots(panels, ncol=length(panels)) +
    plot_annotation(title="Hub TFs at DMR-associated motifs",
      subtitle=paste0("Top ",TOPN," enriched motifs per contrast; STRINGdb PPI; node colour=degree"),
      theme=theme(plot.title=element_text(face="bold", size=13)))
  ggsave(file.path(OUT_DIR,"tf_hub_network_combined.pdf"),
         combined, width=10*length(panels), height=9, device=cairo_pdf)
  message("saved: tf_hub_network_combined.pdf")
}
message("\nDone. Outputs in: ", OUT_DIR)
