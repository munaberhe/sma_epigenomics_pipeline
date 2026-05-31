#!/usr/bin/env Rscript
# 16_tf_motif_enrichment.R
# Transcription factor motif enrichment analysis on ASO-specific DMRs
# Tests whether the 151 ASO-specific DMRs are enriched for binding sites
# of specific TFs compared to matched background regions.
# If ISL1, ROBO2, or NeuroD1 motifs come up enriched this directly supports
# the neural pathway GO enrichment we see in the ASO-specific DMR set.
#
# Method: monaLisa with JASPAR2020 vertebrate TF motifs
# Background: random genomic regions matched for length + chromosome distribution
# Input: results/dmr/dmr_ASO_specific.rds
# Output: results/tf_motif/
#

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TFBSTools)
  library(JASPAR2020)
  library(monaLisa)
  library(ggplot2)
})
.libPaths(c("~/R/library", .libPaths()))

DMR_DIR <- "results/dmr"
OUT_DIR <- "results/tf_motif"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

KEEP_CHRS <- paste0("chr", 1:22)

# Load ASO-specific DMRs -- the 151 DMRs unique to ASO_CTRL contrast
# These are the most biologically interesting set since they represent
# the pure nusinersen off-target epigenetic signature
message("Loading ASO-specific DMRs...")
dmr_specific <- readRDS(file.path(DMR_DIR, "dmr_ASO_specific.rds"))
dmr_specific <- dmr_specific[as.character(seqnames(dmr_specific)) %in% KEEP_CHRS]
message("  ASO-specific DMRs (autosomes): ", length(dmr_specific))

# Make matched background regions -- same chromosome distribution and
# same width distribution as the query DMRs, not overlapping the DMRs
# We use 5x background to give the test more power
make_background <- function(query, n_bg=NULL, seed=42) {
  if (is.null(n_bg)) n_bg <- length(query) * 5
  set.seed(seed)
  chr_sizes <- c(
    chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
    chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
    chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
    chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
    chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
    chr21=46709983, chr22=50818468)
  chr_counts <- table(as.character(seqnames(query)))
  chr_counts <- chr_counts[names(chr_counts) %in% names(chr_sizes)]
  chr_probs  <- chr_counts / sum(chr_counts)
  query_widths <- width(query)
  bg_list <- list()
  attempts <- 0
  while (length(bg_list) < n_bg && attempts < n_bg * 20) {
    attempts <- attempts + 1
    chr <- sample(names(chr_probs), 1, prob=chr_probs)
    w   <- sample(query_widths, 1)
    maxs <- chr_sizes[chr] - w
    if (maxs < 1) next
    s <- sample(1:maxs, 1)
    candidate <- GRanges(chr, IRanges(s, s+w-1))
    if (length(findOverlaps(candidate, query)) == 0)
      bg_list[[length(bg_list)+1]] <- candidate
  }
  message("  Background regions created: ", length(bg_list))
  do.call(c, bg_list)
}

message("Making matched background regions (5x)...")
bg_regions <- make_background(dmr_specific)

# Extract DNA sequences from hg38 for both DMRs and background
# getSeq() pulls the sequence at each GRanges interval from the BSgenome object
message("Extracting DNA sequences from hg38...")
dmr_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, dmr_specific)
bg_seqs  <- getSeq(BSgenome.Hsapiens.UCSC.hg38, bg_regions)
message("  DMR sequences: ", length(dmr_seqs))
message("  Background sequences: ", length(bg_seqs))

# Save FASTA files for reference and potential use with other tools
writeXStringSet(dmr_seqs, file.path(OUT_DIR, "ASO_specific_DMR_sequences.fa"))
writeXStringSet(bg_seqs,  file.path(OUT_DIR, "background_sequences.fa"))
message("  FASTA files saved")

# Load JASPAR2020 vertebrate TF motifs
# Using vertebrate collection which includes human TFs like ISL1, ROBO2 targets
message("Loading JASPAR2020 TF motifs...")
opts     <- list(collection="CORE", tax_group="vertebrates", all_versions=FALSE)
pfm_list <- getMatrixSet(JASPAR2020, opts)
message("  Loaded ", length(pfm_list), " TF motifs from JASPAR2020")
pwm_list <- toPWM(pfm_list)
message("  Converted to PWM format")

# Run motif enrichment using monaLisa calcBinnedMotifEnrR
# Bin 1 = ASO DMR sequences, Bin 0 = background sequences
# Tests which motifs are enriched in DMRs vs background with GC bias correction
message("Running motif enrichment analysis...")
message("  This may take 10-20 minutes...")

all_seqs <- c(dmr_seqs, bg_seqs)
bins <- factor(c(rep("ASO_DMR", length(dmr_seqs)),
                 rep("Background", length(bg_seqs))),
               levels=c("Background", "ASO_DMR"))

se <- calcBinnedMotifEnrR(
  seqs    = all_seqs,
  bins    = bins,
  pwmL    = pwm_list,
  BPPARAM = BiocParallel::MulticoreParam(4)
)

message("  Motif enrichment complete")
saveRDS(se, file.path(OUT_DIR, "motif_enrichment_results.rds"))

# Extract and save results
log_odds    <- assay(se, "log2enr")[, "ASO_DMR"]
padj        <- assay(se, "padj")[, "ASO_DMR"]
motif_names <- rowData(se)$motif.name

results_df <- data.frame(
  motif_id   = rownames(se),
  motif_name = motif_names,
  log2_enr   = log_odds,
  padj       = padj,
  stringsAsFactors = FALSE
)

results_sig <- results_df[!is.na(results_df$padj) & results_df$padj < 0.05, ]
results_sig <- results_sig[order(results_sig$log2_enr, decreasing=TRUE), ]

message("\nSignificant enriched motifs (padj<0.05, log2enr>0): ",
        sum(results_sig$log2_enr > 0))
message("Top 20:")
print(head(results_sig[results_sig$log2_enr > 0, ], 20), row.names=FALSE)

write.csv(results_df[order(results_df$log2_enr, decreasing=TRUE), ],
          file.path(OUT_DIR, "motif_enrichment_all.csv"), row.names=FALSE)
write.csv(results_sig[results_sig$log2_enr > 0, ],
          file.path(OUT_DIR, "motif_enrichment_significant_enriched.csv"),
          row.names=FALSE)

# Check specifically for neural TFs expected from GO enrichment results
neural_tfs <- c("ISL1","ISL2","ROBO","NEUROD1","NEUROD2","ASCL1",
                "NHLH1","NHLH2","OLIG2","NKX2","NKX6","MNX1",
                "HB9","CHAT","SOX10","PAX6","ATOH1")

message("\nChecking for neural TF motifs:")
for (tf in neural_tfs) {
  hits <- results_df[grepl(tf, results_df$motif_name, ignore.case=TRUE), ]
  if (nrow(hits) > 0) {
    cat(sprintf("  %-12s log2enr=%+.3f  padj=%.3e\n",
                tf, hits$log2_enr[1], hits$padj[1]))
  }
}

# Plot top 30 enriched motifs
top_n  <- 30
plot_df <- head(results_sig[results_sig$log2_enr > 0, ], top_n)

if (nrow(plot_df) > 0) {
  plot_df$motif_name <- factor(plot_df$motif_name,
                                levels=rev(plot_df$motif_name))
  plot_df$sig_level  <- cut(-log10(plot_df$padj),
                             breaks=c(0,1.3,2,3,Inf),
                             labels=c("p<0.05","p<0.01","p<0.001","p<0.0001"))
  p <- ggplot(plot_df, aes(x=log2_enr, y=motif_name, fill=sig_level)) +
    geom_col() +
    scale_fill_manual(
      values=c("p<0.05"="#b8d4e8","p<0.01"="#4b9eff",
               "p<0.001"="#1D6FA4","p<0.0001"="#0d3a6e"),
      name="Significance", drop=FALSE) +
    geom_vline(xintercept=0, linewidth=0.5, colour="grey30") +
    labs(title="TF motif enrichment at ASO-specific DMRs",
         subtitle=paste0("151 ASO-specific DMRs vs matched background (5x)\n",
                         "JASPAR2020 vertebrate motifs, monaLisa binned enrichment"),
         x="Log2 enrichment (ASO DMRs vs background)", y=NULL) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"),
          plot.subtitle=element_text(size=9, colour="grey40"),
          axis.text.y=element_text(size=8))
  ggsave(file.path(OUT_DIR, "motif_enrichment_top30.pdf"), p, width=10, height=8)
  message("Saved: motif_enrichment_top30.pdf")
}

# Plot top 30 depleted motifs
top_dep <- head(results_sig[results_sig$log2_enr < 0, ], top_n)
top_dep <- top_dep[order(top_dep$log2_enr), ]

if (nrow(top_dep) > 0) {
  top_dep$motif_name <- factor(top_dep$motif_name, levels=rev(top_dep$motif_name))
  p2 <- ggplot(top_dep, aes(x=log2_enr, y=motif_name)) +
    geom_col(fill="#D94F3D") +
    geom_vline(xintercept=0, linewidth=0.5, colour="grey30") +
    labs(title="TF motifs depleted at ASO-specific DMRs",
         subtitle="Motifs enriched in background vs ASO DMRs",
         x="Log2 enrichment (negative = depleted in DMRs)", y=NULL) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"), axis.text.y=element_text(size=8))
  ggsave(file.path(OUT_DIR, "motif_depletion_top30.pdf"), p2, width=10, height=8)
  message("Saved: motif_depletion_top30.pdf")
}

message("\nDone. Outputs in: ", OUT_DIR)
