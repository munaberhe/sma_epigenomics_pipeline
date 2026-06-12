.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
})
source("scripts/00_sma_palette.R")
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

# ---- H3K9me2 enrichment at ASO DMRs ----
# Marasco et al. 2022 Cell, GSE167762, HEK293T H3K9me2 vs input ChIP-seq.
# Question: do WGBS-called ASO DMRs sit in heterochromatic (H3K9me2) regions
# more than expected by chance?
# Spec from Radu: delta = mean(ASO ChIP) - mean(CTR ChIP) per region, compare
# to a matched set of random regions of identical size (n = 1500).

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
MIN_CYTOSINES <- 6      # drop DMRs with sparse CpG support
N_REGIONS     <- 1500   # match Radu's spec for DMR + background set size
SEED          <- 42

# ---- load DMRs ----
message("loading DMRs...")
dmr_aso      <- readRDS(file.path(DMR_DIR, "dmr_ASO_CTRL_vs_Scramble_CTRL.rds"))
dmr_specific <- readRDS(file.path(DMR_DIR, "dmr_ASO_specific.rds"))
dmr_aso      <- dmr_aso[as.character(seqnames(dmr_aso))      %in% KEEP_CHRS]
dmr_specific <- dmr_specific[as.character(seqnames(dmr_specific)) %in% KEEP_CHRS]

# QC filter: require >= 6 informative cytosines per DMR
if ("cytosinesCount" %in% colnames(mcols(dmr_aso))) {
  before <- length(dmr_aso)
  dmr_aso <- dmr_aso[dmr_aso$cytosinesCount >= MIN_CYTOSINES]
  message("  ASO DMRs cytosine filter: ", before, " -> ", length(dmr_aso))
}
if ("cytosinesCount" %in% colnames(mcols(dmr_specific))) {
  before <- length(dmr_specific)
  dmr_specific <- dmr_specific[dmr_specific$cytosinesCount >= MIN_CYTOSINES]
  message("  ASO-specific cytosine filter: ", before, " -> ", length(dmr_specific))
}
message("  ASO DMRs: ", length(dmr_aso))
message("  ASO-specific DMRs: ", length(dmr_specific))

# ---- helpers ----

# mean bigWig score across each range; NA if range has no coverage
bw_mean_over_ranges <- function(bw_file, gr) {
  message("  importing: ", basename(bw_file))
  bw <- import(bw_file, which=gr, as="NumericList")
  vapply(bw, function(x) if(length(x)==0) NA_real_ else mean(x, na.rm=TRUE),
         numeric(1))
}

# build n_bg random regions matched to query by chromosome and width.
# excludes anything overlapping the query so background is independent.
make_background <- function(query, n_bg=NULL, seed=SEED) {
  if (is.null(n_bg)) n_bg <- length(query)
  set.seed(seed)
  chr_sizes <- c(
    chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
    chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
    chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
    chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
    chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
    chr21=46709983, chr22=50818468)
  chr_counts   <- table(as.character(seqnames(query)))
  chr_counts   <- chr_counts[names(chr_counts) %in% names(chr_sizes)]
  chr_probs    <- chr_counts / sum(chr_counts)
  query_widths <- width(query)
  bg_list  <- list()
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

# ---- subsample DMRs to match background n ----
# Radu: "choose another random set of regions with the same size as dmrs. 1500"
# previous version sampled 5000 ASO DMRs and 1500 BG; n mismatch inflates
# variance and breaks the matched comparison.
set.seed(1)
if (length(dmr_aso) > N_REGIONS) {
  dmr_aso_sub <- dmr_aso[sample(length(dmr_aso), N_REGIONS)]
} else {
  dmr_aso_sub <- dmr_aso
}
if (length(dmr_specific) > N_REGIONS) {
  dmr_specific_sub <- dmr_specific[sample(length(dmr_specific), N_REGIONS)]
} else {
  dmr_specific_sub <- dmr_specific
}
message("subsampled to: ASO=", length(dmr_aso_sub),
        " specific=", length(dmr_specific_sub))

# ---- background regions matched by chr + width ----
message("making background regions...")
bg_aso      <- make_background(dmr_aso_sub,      n_bg=length(dmr_aso_sub))
bg_specific <- make_background(dmr_specific_sub, n_bg=length(dmr_specific_sub))

# ---- extract ChIP signal ----
sig_aso      <- compute_signal(dmr_aso_sub,      "ASO_DMR")
sig_bg_aso   <- compute_signal(bg_aso,           "Background")
sig_specific <- compute_signal(dmr_specific_sub, "ASO_specific_DMR")
sig_bg_spec  <- compute_signal(bg_specific,      "Background")

# mean across replicates per condition
for (df_name in c("sig_aso","sig_bg_aso","sig_specific","sig_bg_spec")) {
  df          <- get(df_name)
  df$mean_CTR <- rowMeans(df[,c("CTR_R1","CTR_R2")], na.rm=TRUE)
  df$mean_ASO <- rowMeans(df[,c("ASO_R1","ASO_R2")], na.rm=TRUE)
  assign(df_name, df)
}

# ---- baseline CTR-only tests (kept for the summary table) ----
test1 <- wilcox.test(sig_aso$mean_CTR,      sig_bg_aso$mean_CTR,  alternative="greater")
test2 <- wilcox.test(sig_specific$mean_CTR, sig_bg_spec$mean_CTR, alternative="greater")
message("Wilcoxon p-values (CTR signal, one-sided):")
message("  ASO DMRs vs BG:     p = ", signif(test1$p.value, 3))
message("  ASO-specific vs BG: p = ", signif(test2$p.value, 3))

summary_df <- data.frame(
  contrast       = c("ASO_DMRs", "ASO_specific_DMRs"),
  n              = c(nrow(sig_aso), nrow(sig_specific)),
  median_DMR_CTR = c(median(sig_aso$mean_CTR,      na.rm=TRUE),
                     median(sig_specific$mean_CTR, na.rm=TRUE)),
  median_BG_CTR  = c(median(sig_bg_aso$mean_CTR,   na.rm=TRUE),
                     median(sig_bg_spec$mean_CTR,  na.rm=TRUE)),
  wilcox_p       = c(signif(test1$p.value,3), signif(test2$p.value,3))
)
print(summary_df, row.names=FALSE)
write.table(summary_df, file.path(OUT_DIR, "h3k9me2_signal_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ---- box panels: raw signal (for reference, not the primary figure) ----
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
               "does ASO raise H3K9me2 at DMR loci?",
               "Mean H3K9me2 signal")

combined <- p1 / p2 +
  plot_annotation(title="H3K9me2 enrichment at WGBS DMR loci",
                  theme=theme(plot.title=element_text(face="bold", size=13)))
ggsave(file.path(OUT_DIR, "h3k9me2_signal_boxplot.pdf"), combined,
       width=10, height=9, device=cairo_pdf)

# ---- primary figure: delta H3K9me2 (ASO - CTR) ----
# DMR delta vs matched-size background delta. Two-sided so we can detect
# either gain or depletion vs background.
message("computing delta H3K9me2 signal (ASO - CTR)...")
delta_dmr <- data.frame(
  delta = sig_aso$mean_ASO - sig_aso$mean_CTR,
  group = "ASO DMRs"
)
delta_bg <- data.frame(
  delta = sig_bg_aso$mean_ASO - sig_bg_aso$mean_CTR,
  group = paste0("Background (n=", nrow(sig_bg_aso), ")")
)
delta_df <- rbind(delta_dmr, delta_bg)
delta_df$group <- factor(delta_df$group,
                         levels=c(unique(delta_bg$group), "ASO DMRs"))

wtest_delta <- wilcox.test(delta_dmr$delta, delta_bg$delta,
                           alternative="two.sided")
message("  delta Wilcoxon (two-sided) p = ", signif(wtest_delta$p.value, 3))
message("  median delta DMR: ", signif(median(delta_dmr$delta, na.rm=TRUE), 3))
message("  median delta BG : ", signif(median(delta_bg$delta,  na.rm=TRUE), 3))

# y-axis clipping: outliers stretch axis to +/- 200 and hide the box bodies.
# 99th percentile of |delta| is ~10-15 so +/- 30 keeps whiskers visible.
y_lim <- c(-15, 15)

p_delta <- ggplot(delta_df, aes(x=group, y=delta, fill=group)) +
  geom_boxplot(outlier.shape=NA, width=0.4) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  scale_fill_manual(values=setNames(c("#cccccc","#1B4F8A"),
                                    levels(delta_df$group))) +
  coord_cartesian(ylim=y_lim) +
  annotate("text", x=1.5, y=y_lim[2]*0.92,
           label=paste0("Wilcoxon two-sided p = ",
                        signif(wtest_delta$p.value, 3)),
           size=3.7, colour="grey20") +
  theme_classic(base_size=12) +
  theme(legend.position="none",
        plot.title=element_text(face="bold"),
        axis.text.x=element_text(size=11)) +
  labs(title="H3K9me2 delta (ASO - CTR): DMRs vs matched background",
       subtitle=paste0("DMRs n=", nrow(delta_dmr),
                       ", background n=", nrow(delta_bg),
                       " (chr/width matched)"),
       x=NULL, y="Delta H3K9me2 signal (ASO - CTR)")

ggsave(file.path(OUT_DIR, "h3k9me2_delta_signal.pdf"),
       p_delta, width=10, height=6, device=cairo_pdf)
message("saved: h3k9me2_delta_signal.pdf")

# ---- SMN2 locus spot check ----
# focal, not genome-wide: report mean signal across replicates per condition.
message("\nSMN2 locus H3K9me2 signal...")
smn2    <- GRanges("chr5", IRanges(70049638, 70078522))
smn2_df <- data.frame(
  condition = names(BW_FILES),
  signal    = sapply(BW_FILES, function(f) bw_mean_over_ranges(f, smn2))
)
print(smn2_df, row.names=FALSE)
write.table(smn2_df, file.path(OUT_DIR, "smn2_h3k9me2_signal.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# summary line for the slide deck
smn2_ctr <- mean(smn2_df$signal[grepl("^CTR", smn2_df$condition)])
smn2_aso <- mean(smn2_df$signal[grepl("^ASO", smn2_df$condition)])
message(sprintf("  SMN2: CTR=%.2f ASO=%.2f fold=%.2fx",
                smn2_ctr, smn2_aso, smn2_aso/smn2_ctr))

message("\ndone. outputs in: ", OUT_DIR)
