#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ChIPseeker)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")


# H3K9me2 signal enrichment at ASO-associated DMRs.
# Data: Marasco et al. 2022 Cell (GSE167762), HEK293T cells.
# Tests whether ASO-linked DMRs overlap with H3K9me2-marked heterochromatin.

BW_FILES <- c(
  CTR_R1 = "data/external/GSE167762_H3K9me2/GSM6063702_CTRvsInp_CTR_R1.bw",
  CTR_R2 = "data/external/GSE167762_H3K9me2/GSM6063706_CTRvsInp_CTR_R2.bw",
  ASO_R1 = "data/external/GSE167762_H3K9me2/GSM6063703_ASOvsInp_ASO_R1.bw",
  ASO_R2 = "data/external/GSE167762_H3K9me2/GSM6063707_ASOvsInp_ASO_R2.bw"
)

DMR_DIR <- "results/dmr"
OUT_DIR <- "results/h3k9me2_overlap"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

KEEP_CHRS <- paste0("chr", 1:22)

# load DMRs
message("loading DMRs...")
dmr_aso      <- readRDS(file.path(DMR_DIR, "dmr_ASO_CTRL_vs_Scramble_CTRL.rds"))
dmr_specific <- readRDS(file.path(DMR_DIR, "dmr_ASO_specific.rds"))
dmr_aso      <- dmr_aso[as.character(seqnames(dmr_aso)) %in% KEEP_CHRS]
dmr_specific <- dmr_specific[as.character(seqnames(dmr_specific)) %in% KEEP_CHRS]
message("  ASO DMRs: ", length(dmr_aso))
message("  ASO-specific DMRs: ", length(dmr_specific))

# mean bigWig signal over a set of ranges
bw_mean_over_ranges <- function(bw_file, gr) {
  message("  importing: ", basename(bw_file))
  bw <- import(bw_file, which=gr, as="NumericList")
  vapply(bw, function(x) if(length(x)==0) NA_real_ else mean(x, na.rm=TRUE),
         numeric(1))
}

# matched background regions with same width/chromosome distribution
make_background <- function(query, n_bg=NULL, seed=42) {
  if (is.null(n_bg)) n_bg <- length(query)
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
  while (length(bg_list) < n_bg && attempts < n_bg*10) {
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
  do.call(c, bg_list)
}

compute_signal <- function(gr, label) {
  message("\ncomputing signal over ", label, " (n=", length(gr), ")...")
  data.frame(
    CTR_R1 = bw_mean_over_ranges(BW_FILES["CTR_R1"], gr),
    CTR_R2 = bw_mean_over_ranges(BW_FILES["CTR_R2"], gr),
    ASO_R1 = bw_mean_over_ranges(BW_FILES["ASO_R1"], gr),
    ASO_R2 = bw_mean_over_ranges(BW_FILES["ASO_R2"], gr),
    group  = label
  )
}

# subsample ASO DMRs for speed
set.seed(1)
dmr_aso_sub <- dmr_aso[sample(length(dmr_aso), min(5000, length(dmr_aso)))]

message("making background regions...")
bg_aso      <- make_background(dmr_aso_sub)
bg_specific <- make_background(dmr_specific, n_bg=length(dmr_specific)*5)

sig_aso      <- compute_signal(dmr_aso_sub,  "ASO_DMR")
sig_bg_aso   <- compute_signal(bg_aso,       "Background")
sig_specific <- compute_signal(dmr_specific, "ASO_specific_DMR")
sig_bg_spec  <- compute_signal(bg_specific,  "Background")

# add mean columns
for (df_name in c("sig_aso","sig_bg_aso","sig_specific","sig_bg_spec")) {
  df          <- get(df_name)
  df$mean_CTR <- rowMeans(df[,c("CTR_R1","CTR_R2")], na.rm=TRUE)
  df$mean_ASO <- rowMeans(df[,c("ASO_R1","ASO_R2")], na.rm=TRUE)
  assign(df_name, df)
}

# Wilcoxon tests
test1 <- wilcox.test(sig_aso$mean_CTR,      sig_bg_aso$mean_CTR,  alternative="greater")
test2 <- wilcox.test(sig_specific$mean_CTR, sig_bg_spec$mean_CTR, alternative="greater")
message("Wilcoxon p-values:")
message("  ASO DMRs vs BG: p = ", signif(test1$p.value, 3))
message("  ASO-specific vs BG: p = ", signif(test2$p.value, 3))

summary_df <- data.frame(
  contrast         = c("ASO_DMRs", "ASO_specific_DMRs"),
  n                = c(nrow(sig_aso), nrow(sig_specific)),
  median_DMR_CTR   = c(median(sig_aso$mean_CTR,      na.rm=TRUE),
                        median(sig_specific$mean_CTR, na.rm=TRUE)),
  median_BG_CTR    = c(median(sig_bg_aso$mean_CTR,   na.rm=TRUE),
                        median(sig_bg_spec$mean_CTR,  na.rm=TRUE)),
  wilcox_p         = c(signif(test1$p.value,3), signif(test2$p.value,3))
)
print(summary_df, row.names=FALSE)
write.table(summary_df, file.path(OUT_DIR, "h3k9me2_signal_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# plots
plot_df <- rbind(
  data.frame(signal=sig_aso$mean_CTR,      group="ASO DMRs",   set="All ASO DMRs"),
  data.frame(signal=sig_bg_aso$mean_CTR,   group="Background", set="All ASO DMRs"),
  data.frame(signal=sig_specific$mean_CTR, group="ASO DMRs",   set="ASO-specific DMRs"),
  data.frame(signal=sig_bg_spec$mean_CTR,  group="Background", set="ASO-specific DMRs")
)
saveRDS(plot_df, file.path(OUT_DIR, "h3k9me2_plot_df.rds"))
plot_df$group <- factor(plot_df$group, levels=c("Background","ASO DMRs"))

plot_df2 <- rbind(
  data.frame(signal=sig_aso$mean_CTR,      condition="CTR", set="All ASO DMRs"),
  data.frame(signal=sig_aso$mean_ASO,      condition="ASO", set="All ASO DMRs"),
  data.frame(signal=sig_specific$mean_CTR, condition="CTR", set="ASO-specific DMRs"),
  data.frame(signal=sig_specific$mean_ASO, condition="ASO", set="ASO-specific DMRs")
)
plot_df2$condition <- factor(plot_df2$condition, levels=c("CTR","ASO"))

make_box <- function(df, x, fill_col, title, subtitle, y_lab) {
  ggplot(df, aes_string(x=x, y="signal", fill=x)) +
    geom_boxplot(outlier.size=0.3, outlier.alpha=0.3) +
    facet_wrap(~set, scales="free_y") +
    scale_fill_manual(values=fill_col) +
    theme_classic(base_size=11) +
    theme(legend.position="none", strip.text=element_text(face="bold")) +
    labs(title=title, subtitle=subtitle, x=NULL, y=y_lab)
}

p1 <- make_box(plot_df, "group",
               c("Background"="#cccccc","ASO DMRs"="#1B4F8A"),
               "H3K9me2 at ASO DMRs vs background",
               "CTR ChIP-seq (Marasco et al. 2022)", "Mean H3K9me2 signal")
p2 <- make_box(plot_df2, "condition",
               c("CTR"="#6B7280","ASO"="#1B4F8A"),
               "H3K9me2 at DMRs: CTR vs ASO",
               "tests whether ASO increases H3K9me2 at DMR loci",
               "Mean H3K9me2 signal")

combined <- p1 / p2 +
  plot_annotation(title="H3K9me2 enrichment at WGBS DMR loci",
                  theme=theme(plot.title=element_text(face="bold", size=13)))
ggsave(file.path(OUT_DIR, "h3k9me2_signal_boxplot.pdf"), combined, width=10, height=9)

# SMN2 locus spot check
message("\nSMN2 locus H3K9me2 signal...")
smn2    <- GRanges("chr5", IRanges(70049638, 70078522))
smn2_df <- data.frame(
  condition = names(BW_FILES),
  signal    = sapply(BW_FILES, function(f) bw_mean_over_ranges(f, smn2))
)
print(smn2_df, row.names=FALSE)
write.table(smn2_df, file.path(OUT_DIR, "smn2_h3k9me2_signal.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("\ndone. outputs in: ", OUT_DIR)

OUT <- "results/smn2_enhancer"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# SMN2 gene body and extended window
SMN2_START  <- 70049638
SMN2_END    <- 70078522
WINDOW_START <- 69950000
WINDOW_END   <- 70150000
CHR <- "chr5"

# Exon boundaries (Alberto's convention)
EXONS <- data.frame(
  exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
  start = c(70049638,70053107,70056229,70063044,70069090,
            70069235,70070641,70076521,70077019),
  end   = c(70050437,70053264,70056357,70063153,70069186,
            70069330,70070751,70076574,70077592)
)

# Load narrowPeak files
load_peaks <- function(file, label) {
  df <- read.table(gzfile(file), header=FALSE, sep="\t",
    col.names=c("chr","start","end","name","score",
                "strand","fc","pval","qval","summit"))
  df <- df[df$chr==CHR & df$start>=WINDOW_START & df$end<=WINDOW_END,]
  df$label <- label
  df$mid <- (df$start + df$end) / 2
  df
}

peaks_ctrl1 <- load_peaks(
  "data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep1.narrowPeak.gz",
  "CTRL Rep1")
peaks_ctrl2 <- load_peaks(
  "data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep2.narrowPeak.gz",
  "CTRL Rep2")
peaks_vpa1  <- load_peaks(
  "data/external/h3k27ac_gse246399/H3K27ac_VPA_Rep1.narrowPeak.gz",
  "VPA Rep1")

all_peaks <- rbind(peaks_ctrl1, peaks_ctrl2, peaks_vpa1)
all_peaks$label <- factor(all_peaks$label,
  levels=c("CTRL Rep1","CTRL Rep2","VPA Rep1"))

message("Peaks found:")
print(table(all_peaks$label))

# Save peak table
write.csv(all_peaks[,c("chr","start","end","fc","qval","label")],
  file.path(OUT, "smn2_h3k27ac_peaks.csv"), row.names=FALSE)
message("Saved: smn2_h3k27ac_peaks.csv")

# ── FIGURE: H3K27ac peaks at SMN2 locus ──────────────────────
COLS <- c("CTRL Rep1"="#1D6FA4", "CTRL Rep2"="#2E9B6F", "VPA Rep1"="#F0A500")

# Panel A: peak map showing position and strength
p_peaks <- ggplot(all_peaks,
    aes(xmin=start/1e6, xmax=end/1e6,
        ymin=as.numeric(label)-0.35,
        ymax=as.numeric(label)+0.35,
        fill=label)) +
  geom_rect(alpha=0.8) +
  # Gene body
  annotate("rect", xmin=SMN2_START/1e6, xmax=SMN2_END/1e6,
           ymin=0.3, ymax=0.7, fill="grey80", colour="grey50", linewidth=0.3) +
  annotate("text", x=(SMN2_START+SMN2_END)/2/1e6, y=0.5,
           label="SMN2", size=3, fontface="bold") +
  annotate("text", x=70.025, y=3.6,
           label="SERF1B", size=2.5, colour="grey40", fontface="italic") +
  # Exons
  geom_rect(data=EXONS,
    aes(xmin=start/1e6, xmax=end/1e6, ymin=0.2, ymax=0.8),
    inherit.aes=FALSE, fill="#1B4F8A", colour="white", linewidth=0.2) +
  # E7 in red
  annotate("rect",
    xmin=70076521/1e6, xmax=70076574/1e6,
    ymin=0.15, ymax=0.85, fill="#D94F3D", colour="white", linewidth=0.2) +
  annotate("text", x=70076547/1e6, y=0.05,
           label="E7\n(ASO target)", size=2, colour="#D94F3D") +
  scale_fill_manual(values=COLS, name=NULL) +
  scale_y_continuous(breaks=1:3,
    labels=c("CTRL Rep1","CTRL Rep2","VPA Rep1")) +
  scale_x_continuous(labels=function(x) paste0(x, " Mb")) +
  coord_cartesian(xlim=c(70.0, 70.085)) +
  labs(title="(A) H3K27ac ChIP-seq peaks at SMN2 locus",
       subtitle=NULL,
       x="chr5 position", y=NULL) +
  theme_bw(base_size=10) +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=8),
        panel.grid.minor=element_blank())

# Panel B: fold enrichment at each peak coloured by condition
p_fc <- ggplot(all_peaks,
    aes(x=mid/1e6, y=fc, colour=label, shape=label)) +
  geom_point(size=3, alpha=0.9) +
  geom_vline(xintercept=c(SMN2_START/1e6, SMN2_END/1e6),
             linetype="dashed", colour="grey60", linewidth=0.4) +
  annotate("text", x=SMN2_START/1e6, y=max(all_peaks$fc)*0.95,
           label="SMN2\n5'", size=2.5, hjust=1.1, colour="grey40") +
  annotate("text", x=SMN2_END/1e6, y=max(all_peaks$fc)*0.95,
           label="SMN2\n3'", size=2.5, hjust=-0.1, colour="grey40") +
  annotate("rect", xmin=70076521/1e6, xmax=70076574/1e6,
           ymin=-Inf, ymax=Inf, fill="#D94F3D", alpha=0.1) +
  scale_colour_manual(values=COLS, name=NULL) +
  scale_shape_manual(values=c(16,17,15), name=NULL) +
  scale_x_continuous(labels=function(x) paste0(x, " Mb")) +
  coord_cartesian(xlim=c(70.0, 70.085)) +
  labs(title="(B) H3K27ac fold enrichment by peak position",
       subtitle=NULL,
       x="chr5 position", y="Fold enrichment over input") +
  theme_bw(base_size=10) +
  theme(legend.position="top",
        plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=8),
        panel.grid.minor=element_blank())

fig <- (p_peaks / p_fc) +
  plot_annotation(
    title="SMN2 locus H3K27ac enhancer analysis",
    subtitle=NULL,
    theme=theme(plot.title=element_text(face="bold", size=11),
                plot.subtitle=element_text(size=8, colour="grey30"))
  )

ggsave(file.path(OUT, "SMN2_H3K27ac_enhancer_analysis.pdf"),
       fig, width=14, height=10)
message("Saved: SMN2_H3K27ac_enhancer_analysis")

# Summary table for thesis
summary_df <- all_peaks[,c("chr","start","end","fc","qval","label")]
summary_df$position <- ifelse(summary_df$end < SMN2_START, "upstream",
                       ifelse(summary_df$start > SMN2_END, "downstream",
                       ifelse(summary_df$start < SMN2_START+5000, "promoter/5-end",
                       "gene body")))
write.csv(summary_df, file.path(OUT, "smn2_h3k27ac_peak_summary.csv"),
  row.names=FALSE)
message("\nPeak summary:")
print(table(summary_df$label, summary_df$position))
message("\nDone. Results in: ", OUT)


# ── ITEM 9: SMN2 introns 6 and 7 H3K27ac check ───────────────────────────────
message('=== ITEM 9: SMN2 intron 6-7 H3K27ac analysis ===')

# SMN2 exon coordinates (Alberto's convention, E7 = penultimate exon)
# Intron 6 = between E6 and E7
# Intron 7 = between E7 and E8
EXONS <- data.frame(
  exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
  start = c(70049638,70053107,70056229,70063044,70069090,
            70069235,70070641,70076521,70077019),
  end   = c(70050437,70053264,70056357,70063153,70069186,
            70069330,70070751,70076574,70077592)
)

# Intron boundaries
intron6_start <- EXONS$end[EXONS$exon=="E6"]   # after E6
intron6_end   <- EXONS$start[EXONS$exon=="E7"] # before E7
intron7_start <- EXONS$end[EXONS$exon=="E7"]   # after E7
intron7_end   <- EXONS$start[EXONS$exon=="E8"] # before E8

message("Intron 6: chr5:", intron6_start, "-", intron6_end,
        " (", intron6_end-intron6_start, "bp)")
message("Intron 7: chr5:", intron7_start, "-", intron7_end,
        " (", intron7_end-intron7_start, "bp)")

# Load H3K27ac peaks
load_peaks <- function(file, label) {
  df <- read.table(gzfile(file), header=FALSE, sep='\t',
    col.names=c('chr','start','end','name','score',
                'strand','fc','pval','qval','summit'))
  df$label <- label
  df
}

peaks_ctrl1 <- load_peaks(
  'data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep1.narrowPeak.gz', 'CTRL_Rep1')
peaks_ctrl2 <- load_peaks(
  'data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep2.narrowPeak.gz', 'CTRL_Rep2')
peaks_vpa1  <- load_peaks(
  'data/external/h3k27ac_gse246399/H3K27ac_VPA_Rep1.narrowPeak.gz',  'VPA_Rep1')
all_peaks <- rbind(peaks_ctrl1, peaks_ctrl2, peaks_vpa1)

# Check intron 6
message('\n--- Intron 6 (E6-E7 boundary, flanks ASO target) ---')
int6 <- all_peaks[all_peaks$chr=='chr5' &
                  all_peaks$start >= intron6_start &
                  all_peaks$end   <= intron6_end, ]
if (nrow(int6) > 0) {
  message('PEAKS FOUND in intron 6:')
  print(int6[,c('chr','start','end','fc','qval','label')])
} else {
  message('No H3K27ac peaks in intron 6 in any condition')
}

# Check intron 7
message('\n--- Intron 7 (E7-E8 boundary, 3 prime end) ---')
int7 <- all_peaks[all_peaks$chr=='chr5' &
                  all_peaks$start >= intron7_start &
                  all_peaks$end   <= intron7_end, ]
if (nrow(int7) > 0) {
  message('PEAKS FOUND in intron 7:')
  print(int7[,c('chr','start','end','fc','qval','label')])
} else {
  message('No H3K27ac peaks in intron 7 in any condition')
}

# Extended check — 5kb around E7 (ASO target region)
message('\n--- 5kb window around E7 (ASO target ±5kb) ---')
e7_window <- all_peaks[all_peaks$chr=='chr5' &
                        all_peaks$start >= (EXONS$start[EXONS$exon=="E7"] - 5000) &
                        all_peaks$end   <= (EXONS$end[EXONS$exon=="E7"]   + 5000), ]
if (nrow(e7_window) > 0) {
  message('PEAKS near E7:')
  print(e7_window[,c('chr','start','end','fc','qval','label')])
} else {
  message('No H3K27ac peaks within 5kb of E7')
}

# Save intron results
intron_results <- data.frame(
  region = c('Intron 6', 'Intron 7', 'E7 ±5kb'),
  coordinates = c(
    paste0('chr5:', intron6_start, '-', intron6_end),
    paste0('chr5:', intron7_start, '-', intron7_end),
    paste0('chr5:', EXONS$start[EXONS$exon=="E7"]-5000, '-',
                    EXONS$end[EXONS$exon=="E7"]+5000)
  ),
  peaks_CTRL_Rep1 = c(
    sum(int6$label=='CTRL_Rep1'),
    sum(int7$label=='CTRL_Rep1'),
    sum(e7_window$label=='CTRL_Rep1')
  ),
  peaks_CTRL_Rep2 = c(
    sum(int6$label=='CTRL_Rep2'),
    sum(int7$label=='CTRL_Rep2'),
    sum(e7_window$label=='CTRL_Rep2')
  ),
  peaks_VPA_Rep1 = c(
    sum(int6$label=='VPA_Rep1'),
    sum(int7$label=='VPA_Rep1'),
    sum(e7_window$label=='VPA_Rep1')
  )
)
write.csv(intron_results,
  'results/smn2_enhancer/SMN2_intron67_H3K27ac_summary.csv',
  row.names=FALSE)
message('\nSaved: SMN2_intron67_H3K27ac_summary.csv')

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

