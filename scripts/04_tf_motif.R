#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TFBSTools)
  library(JASPAR2020)
  library(monaLisa)
  library(ggplot2)
})

# TF motif enrichment at ASO-specific DMRs.
# Uses monaLisa with JASPAR2020 vertebrate motifs.
# Background: random genomic regions matched for length and chromosome distribution.
# If ISL1, ROBO2, NeuroD1 come up enriched this supports the neural GO enrichment.

DMR_DIR <- "results/dmr"
OUT_DIR <- "results/tf_motif"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

KEEP_CHRS <- paste0("chr", 1:22)

message("loading ASO-specific DMRs...")
dmr_specific <- readRDS(file.path(DMR_DIR, "dmr_ASO_specific.rds"))
dmr_specific <- dmr_specific[as.character(seqnames(dmr_specific)) %in% KEEP_CHRS]
message("  n = ", length(dmr_specific))

# matched background regions — same chr distribution, same widths, no overlap with DMRs
make_background <- function(query, n_bg=NULL, seed=42) {
  if (is.null(n_bg)) n_bg <- length(query) * 5
  set.seed(seed)
  chr_sizes <- c(
    chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
    chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
    chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
    chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
    chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
    chr21=46709983,  chr22=50818468)
  chr_counts <- table(as.character(seqnames(query)))
  chr_counts <- chr_counts[names(chr_counts) %in% names(chr_sizes)]
  chr_probs  <- chr_counts / sum(chr_counts)
  query_widths <- width(query)
  bg_list <- list()
  attempts <- 0
  while (length(bg_list) < n_bg && attempts < n_bg*20) {
    attempts <- attempts + 1
    chr  <- sample(names(chr_probs), 1, prob=chr_probs)
    w    <- sample(query_widths, 1)
    maxs <- chr_sizes[chr] - w
    if (maxs < 1) next
    s <- sample(1:maxs, 1)
    candidate <- GRanges(chr, IRanges(s, s+w-1))
    if (length(findOverlaps(candidate, query)) == 0)
      bg_list[[length(bg_list)+1]] <- candidate
  }
  message("  background regions: ", length(bg_list))
  do.call(c, bg_list)
}

message("making background regions (5x)...")
bg_regions <- make_background(dmr_specific)

# extract sequences from hg38
message("extracting sequences...")
dmr_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, dmr_specific)
bg_seqs  <- getSeq(BSgenome.Hsapiens.UCSC.hg38, bg_regions)
writeXStringSet(dmr_seqs, file.path(OUT_DIR, "ASO_specific_DMR_sequences.fa"))
writeXStringSet(bg_seqs,  file.path(OUT_DIR, "background_sequences.fa"))

# load JASPAR2020 vertebrate motifs
message("loading JASPAR2020 motifs...")
pfm_list <- getMatrixSet(JASPAR2020,
  list(collection="CORE", tax_group="vertebrates", all_versions=FALSE))
pwm_list <- toPWM(pfm_list)
message("  ", length(pfm_list), " motifs loaded")

# run motif enrichment — this takes 10-20 min
message("running motif enrichment (10-20 min)...")
all_seqs <- c(dmr_seqs, bg_seqs)
bins <- factor(c(rep("ASO_DMR",    length(dmr_seqs)),
                 rep("Background", length(bg_seqs))),
               levels=c("Background","ASO_DMR"))
se <- calcBinnedMotifEnrR(
  seqs    = all_seqs,
  bins    = bins,
  pwmL    = pwm_list,
  BPPARAM = BiocParallel::MulticoreParam(4)
)
saveRDS(se, file.path(OUT_DIR, "motif_enrichment_results.rds"))

# extract results
results_df <- data.frame(
  motif_id   = rownames(se),
  motif_name = rowData(se)$motif.name,
  log2_enr   = assay(se, "log2enr")[, "ASO_DMR"],
  padj       = assay(se, "padj")[,   "ASO_DMR"],
  stringsAsFactors = FALSE
)
results_sig <- results_df[!is.na(results_df$padj) & results_df$padj < 0.05, ]
results_sig <- results_sig[order(results_sig$log2_enr, decreasing=TRUE), ]

message("\nsignificant enriched motifs: ", sum(results_sig$log2_enr > 0))
print(head(results_sig[results_sig$log2_enr > 0, ], 20), row.names=FALSE)

write.csv(results_df[order(results_df$log2_enr, decreasing=TRUE), ],
          file.path(OUT_DIR, "motif_enrichment_all.csv"), row.names=FALSE)
write.csv(results_sig[results_sig$log2_enr > 0, ],
          file.path(OUT_DIR, "motif_enrichment_significant_enriched.csv"),
          row.names=FALSE)

# spot check for neural TFs expected from GO enrichment
neural_tfs <- c("ISL1","ISL2","ROBO","NEUROD1","NEUROD2","ASCL1",
                "NHLH1","NHLH2","OLIG2","NKX2","NKX6","MNX1",
                "HB9","CHAT","SOX10","PAX6","ATOH1")
message("\nneural TF motifs:")
for (tf in neural_tfs) {
  hits <- results_df[grepl(tf, results_df$motif_name, ignore.case=TRUE), ]
  if (nrow(hits) > 0)
    cat(sprintf("  %-12s log2enr=%+.3f  padj=%.3e\n",
                tf, hits$log2_enr[1], hits$padj[1]))
}

# plot top enriched
top_n   <- 30
plot_df <- head(results_sig[results_sig$log2_enr > 0, ], top_n)
if (nrow(plot_df) > 0) {
  plot_df$motif_name <- factor(plot_df$motif_name, levels=rev(plot_df$motif_name))
  plot_df$sig_level  <- cut(-log10(plot_df$padj),
                             breaks=c(0,1.3,2,3,Inf),
                             labels=c("p<0.05","p<0.01","p<0.001","p<0.0001"))
  p <- ggplot(plot_df, aes(x=log2_enr, y=motif_name, fill=sig_level)) +
    geom_col() +
    scale_fill_manual(values=c("p<0.05"="#b8d4e8","p<0.01"="#4b9eff",
                               "p<0.001"="#1D6FA4","p<0.0001"="#0d3a6e"),
                      name="significance", drop=FALSE) +
    geom_vline(xintercept=0, linewidth=0.5, colour="grey30") +
    labs(title="TF motif enrichment at ASO-specific DMRs",
         subtitle="JASPAR2020 vertebrate motifs, monaLisa binned enrichment",
         x="log2 enrichment (ASO DMRs vs background)", y=NULL) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"), axis.text.y=element_text(size=8))
  ggsave(file.path(OUT_DIR, "motif_enrichment_top30.pdf"), p, width=10, height=8)
  message("saved enrichment plot")
}

# plot top depleted
top_dep <- head(results_sig[results_sig$log2_enr < 0, ], top_n)
top_dep <- top_dep[order(top_dep$log2_enr), ]
if (nrow(top_dep) > 0) {
  top_dep$motif_name <- factor(top_dep$motif_name, levels=rev(top_dep$motif_name))
  p2 <- ggplot(top_dep, aes(x=log2_enr, y=motif_name)) +
    geom_col(fill="#D94F3D") +
    geom_vline(xintercept=0, linewidth=0.5, colour="grey30") +
    labs(title="TF motifs depleted at ASO-specific DMRs",
         x="log2 enrichment (negative = depleted)", y=NULL) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"), axis.text.y=element_text(size=8))
  ggsave(file.path(OUT_DIR, "motif_depletion_top30.pdf"), p2, width=10, height=8)
  message("saved depletion plot")
}


log_odds    <- as.numeric(assay(se, "log2enr")[, "ASO_DMR"])
neg_log_p   <- as.numeric(assay(se, "negLog10Padj")[, "ASO_DMR"])
padj        <- 10^(-neg_log_p)
motif_names <- as.character(rowData(se)[["motif.name"]])

results_df <- data.frame(
  motif_id     = rownames(se),
  motif_name   = motif_names,
  log2_enr     = log_odds,
  negLog10padj = neg_log_p,
  padj         = padj,
  stringsAsFactors = FALSE
)
results_df <- results_df[order(results_df$log2_enr, decreasing=TRUE), ]

message("\ntop 20 enriched motifs:")
print(head(results_df[, c("motif_name","log2_enr","padj")], 20), row.names=FALSE)
message("min padj: ", round(min(results_df$padj, na.rm=TRUE), 4))
message("motifs padj<0.05: ", sum(results_df$padj < 0.05, na.rm=TRUE))

write.csv(results_df, file.path(OUT_DIR, "motif_enrichment_all.csv"), row.names=FALSE)

# spot check neural TFs
neural_tfs <- c("ISL1","ISL2","ROBO","NEUROD1","NEUROD2","ASCL1",
                "NHLH1","NHLH2","OLIG2","NKX2","NKX6","MNX1",
                "HB9","CHAT","SOX10","PAX6","ATOH1")
message("\nneural TF motifs:")
for (tf in neural_tfs) {
  hits <- results_df[grepl(tf, results_df$motif_name, ignore.case=TRUE), ]
  if (nrow(hits) > 0)
    cat(sprintf("  %-12s log2enr=%+.3f  padj=%.3e\n",
                tf, hits$log2_enr[1], hits$padj[1]))
}

# top 20 nominal hits — no significance filter, honest about negative result
top_df <- head(results_df, 20)
top_df$motif_name <- factor(top_df$motif_name, levels=rev(top_df$motif_name))

p <- ggplot(top_df, aes(x=log2_enr, y=motif_name)) +
  geom_col(fill="#4b9eff", alpha=0.8) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey30") +
  labs(title="TF motif enrichment at ASO-specific DMRs (top 20 nominal)",
       x="log2 enrichment (ASO DMRs vs background)", y=NULL,
) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"),
        axis.text.y=element_text(size=9))
ggsave(file.path(OUT_DIR, "motif_enrichment_top20_nominal.pdf"), p, width=10, height=7)
message("saved: motif_enrichment_top20_nominal.pdf")

# volcano plot
results_df$highlight <- ifelse(
  grepl("ISL1|ISL2|NEUROD|ROBO|ASCL|OLIG|NKX|MNX|SOX10",
        results_df$motif_name, ignore.case=TRUE), "neural", "other")

p2 <- ggplot(results_df, aes(x=log2_enr, y=negLog10padj)) +
  geom_point(aes(colour=highlight), size=1.5, alpha=0.7) +
  scale_colour_manual(values=c(neural="#2E9B6F", other="#cccccc"),
                      labels=c(neural="neural TFs", other="other")) +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.20))) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed",
             colour="grey50", linewidth=0.5) +
  annotate("text", x=max(results_df$log2_enr)*0.8,
           y=-log10(0.05)+0.05, label="padj=0.05", size=3, colour="grey50") +
  geom_text(data=head(results_df, 5),
            aes(label=motif_name), size=2.8, hjust=-0.1, colour="#1D6FA4") +
  labs(title="TF motif enrichment — volcano",
       x="log2 enrichment", y="-log10 adjusted p-value") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"))
ggsave(file.path(OUT_DIR, "motif_enrichment_volcano.pdf"), p2, width=10, height=7)
message("saved: motif_enrichment_volcano.pdf")
message("\ndone. outputs in: ", OUT_DIR)
