.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
  library(ggraph)
  library(ggplot2)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/synergy_scan"

SYNERGY_GENES <- c("RELL2", "RNA5S13", "DDIT4L", "MRPS2", "TCEAL4", "KIAA1656", "GNG14")

message("Initialising STRINGdb...")
string_db <- STRINGdb$new(version="12.0", species=9606,
                          score_threshold=200, input_directory="")

message("Mapping ", length(SYNERGY_GENES), " genes to STRING IDs...")
gene_df <- data.frame(gene = SYNERGY_GENES)
mapped <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
message("  mapped: ", nrow(mapped), " of ", length(SYNERGY_GENES))
print(mapped)

if (nrow(mapped) >= 2) {
  message("\nQuerying STRING interactions among these genes...")
  interactions <- string_db$get_interactions(mapped$STRING_id)
  message("  direct interactions found: ", nrow(interactions))
  if (nrow(interactions) > 0) {
    print(interactions)
  } else {
    message("  no direct STRING interactions among these 7 genes at score threshold 200")
  }

  # also check at a much lower confidence threshold, since absence at 200
  # could just mean weak/no evidence, not necessarily nothing at all
  string_db_low <- STRINGdb$new(version="12.0", species=9606,
                                score_threshold=150, input_directory="")
  interactions_low <- string_db_low$get_interactions(mapped$STRING_id)
  message("  interactions at lower threshold (150): ", nrow(interactions_low))
  if (nrow(interactions_low) > 0) print(interactions_low)
} else {
  message("Fewer than 2 genes mapped -- cannot test for a network among them")
}

# plot whatever network exists, even if empty, to confirm visually
png(file.path(OUT_DIR, "synergy_genes_string_network.png"), width=1200, height=900, res=150)
tryCatch({
  string_db$plot_network(mapped$STRING_id)
}, error = function(e) {
  message("plot_network failed (likely because there are no edges): ", e$message)
  plot.new()
  text(0.5, 0.5, "No STRING network edges found among\nthe 7 validated synergy genes", cex=1.3)
})
dev.off()
message("\nsaved: synergy_genes_string_network.png")

message("\nDone.")
