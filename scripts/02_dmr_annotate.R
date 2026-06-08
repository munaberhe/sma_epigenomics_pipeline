#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggplot2)
})

# Annotate DMRs with genomic features and run GO/KEGG enrichment.
# Runs on the genome-wide RDS files produced by 06b_dmrcaller_combine_chr.R.
# Enrichment is split by direction (hypo/hyper) to avoid mixing signals.

OUT_DIR <- "results/dmr_annotation"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

CONTRASTS <- c(
  "ASO_VPA_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_CTRL_vs_Scramble_CTRL",
  "ASO_VPA_vs_ASO_CTRL",
  "ASO_VPA_vs_Scramble_VPA"
)

for (contrast in CONTRASTS) {
  message("\nannotating: ", contrast)
  rds_path <- file.path("results/dmr", paste0("dmr_", contrast, ".rds"))
  if (!file.exists(rds_path)) { message("missing: ", rds_path); next }

  dmrs <- readRDS(rds_path)
  message("DMRs loaded: ", length(dmrs))

  # match chromosome naming to TxDb (UCSC: chr1, chr2...)
  GenomeInfoDb::seqlevelsStyle(dmrs) <- "UCSC"

  # high-confidence filter: at least 6 CpGs per DMR
  dmrs_hc <- dmrs[dmrs$cytosinesCount >= 6]
  message("high-confidence (>=6 CpGs): ", length(dmrs_hc),
          " | hypo: ", sum(dmrs_hc$regionType == "gain"),
          " | hyper: ", sum(dmrs_hc$regionType == "loss"))

  # annotate all DMRs combined
  anno    <- annotatePeak(dmrs_hc, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Hs.eg.db")
  anno_df <- as.data.frame(anno)
  write.csv(anno_df,
            file.path(OUT_DIR, paste0(contrast, "_annotated.csv")),
            row.names=FALSE)

  # top 10 most significant DMRs
  top10 <- anno_df[order(anno_df$pValue), ]
  top10 <- top10[!is.na(top10$SYMBOL), ]
  top10 <- head(top10[, c("seqnames","start","end","cytosinesCount",
                           "regionType","proportion1","proportion2",
                           "pValue","annotation","SYMBOL","GENENAME")], 10)
  top10$direction         <- ifelse(top10$regionType == "gain", "hypo", "hyper")
  top10$methylation_change <- round(top10$proportion2 - top10$proportion1, 3)
  write.csv(top10,
            file.path(OUT_DIR, paste0(contrast, "_top10_genes.csv")),
            row.names=FALSE)
  print(top10[, c("SYMBOL","direction","pValue","methylation_change","annotation")])

  # annotation plots
  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_pie.pdf")), width=8, height=6)
  plotAnnoPie(anno, main=paste0("Genomic features\n", contrast))
  dev.off()

  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_bar.pdf")), width=10, height=6)
  plotAnnoBar(anno, title=paste0("Genomic features — ", contrast))
  dev.off()

  pdf(file.path(OUT_DIR, paste0(contrast, "_TSS_distance.pdf")), width=10, height=6)
  plotDistToTSS(anno, title=paste0("Distance to TSS — ", contrast))
  dev.off()

  # GO and KEGG enrichment — run separately for hypo and hyper
  for (direction in c("gain", "loss")) {
    label <- if (direction == "gain") "hypo" else "hyper"
    sub   <- dmrs_hc[dmrs_hc$regionType == direction]
    message("  ", label, ": ", length(sub), " DMRs")
    if (length(sub) < 10) { message("  too few, skipping"); next }

    anno_sub <- annotatePeak(sub, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
    genes    <- unique(as.data.frame(anno_sub)$geneId)
    genes    <- genes[!is.na(genes)]
    message("  genes: ", length(genes))
    if (length(genes) < 3) { message("  too few genes"); next }

    # GO biological process
    go <- enrichGO(gene=genes, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                   ont="BP", pAdjustMethod="BH",
                   pvalueCutoff=0.05, qvalueCutoff=0.2)
    if (!is.null(go) && nrow(go) > 0) {
      message("  top GO BP: ", go$Description[1])
      write.csv(as.data.frame(go),
                file.path(OUT_DIR, paste0(contrast, "_GO_BP_", label, ".csv")),
                row.names=FALSE)
      pdf(file.path(OUT_DIR, paste0(contrast, "_GO_BP_", label, "_dotplot.pdf")),
          width=10, height=8)
      print(dotplot(go, showCategory=15,
                    title=paste0("GO BP (", label, ") — ", contrast)))
      dev.off()
    } else {
      message("  no significant GO BP terms")
    }

    # KEGG pathways
    kegg <- enrichKEGG(gene=genes, organism="hsa",
                       pvalueCutoff=0.1, qvalueCutoff=0.3, minGSSize=5)
    if (!is.null(kegg) && nrow(kegg) > 0) {
      message("  top KEGG: ", kegg$Description[1])
      write.csv(as.data.frame(kegg),
                file.path(OUT_DIR, paste0(contrast, "_KEGG_", label, ".csv")),
                row.names=FALSE)
      pdf(file.path(OUT_DIR, paste0(contrast, "_KEGG_", label, "_dotplot.pdf")),
          width=10, height=8)
      print(dotplot(kegg, showCategory=15,
                    title=paste0("KEGG (", label, ") — ", contrast)))
      dev.off()
    } else {
      message("  no significant KEGG")
    }
  }
}


CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond_a="ASO_CTRL", cond_b="Scramble_CTRL",
       label="ASO alone vs Scramble CTRL",
       colour=c(ASO_CTRL="#1B4F8A", Scramble_CTRL="#6B7280")),
  list(name="ASO_VPA_vs_Scramble_CTRL",
       cond_a="ASO_VPA", cond_b="Scramble_CTRL",
       label="ASO+VPA vs Scramble CTRL",
       colour=c(ASO_VPA="#B2182B", Scramble_CTRL="#6B7280")),
  list(name="ASO_VPA_vs_ASO_CTRL",
       cond_a="ASO_VPA", cond_b="ASO_CTRL",
       label="VPA effect on ASO background",
       colour=c(ASO_VPA="#B2182B", ASO_CTRL="#1B4F8A")),
  list(name="ASO_VPA_vs_Scramble_VPA",
       cond_a="ASO_VPA", cond_b="Scramble_VPA",
       label="ASO effect on VPA background",
       colour=c(ASO_VPA="#B2182B", Scramble_VPA="#F0A500")),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       cond_a="Scramble_VPA", cond_b="Scramble_CTRL",
       label="VPA alone vs Scramble CTRL",
       colour=c(Scramble_VPA="#F0A500", Scramble_CTRL="#6B7280"))
)

CONDITIONS <- list(
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = paste0("ASO_VPA_",       1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)

KEY_LOCI <- list(
  list(name="SMN2", chr="chr5",
       start=70040000, end=70090000,
       label="SMN2 locus (masked alignment)"),
  list(name="RNA45SN2_promoter", chr="chr21",
       start=8158909, end=8259208,
       label="RNA45SN2 promoter (top hit)"),
  list(name="MTA1-DT", chr="chr14",
       start=105365123, end=105465422,
       label="MTA1-DT locus"),
  list(name="MYO1D", chr="chr17",
       start=32772454, end=32872753,
       label="MYO1D locus")
)

# load pooled methylation data — cached to RDS after first run (~20-30 min)
meth_cache <- file.path(DMR_DIR, "meth_pooled_cache.rds")
if (file.exists(meth_cache)) {
  message("loading cached pooled data...")
  meth_pooled <- readRDS(meth_cache)
} else {
  message("pooling from CpG reports (takes ~20-30 min)...")
  meth_pooled <- list()
  for (cond in names(CONDITIONS)) {
    message("  pooling: ", cond)
    glist <- GRangesList(lapply(CONDITIONS[[cond]], function(s) {
      chr_files <- file.path(BY_CHR_DIR,
        paste0(s, "_", CHROMS, ".CpG_report.txt.gz"))
      chr_files <- chr_files[file.exists(chr_files)]
      do.call(c, lapply(chr_files, readBismark))
    }))
    meth_pooled[[cond]] <- poolMethylationDatasets(glist)
    message("    CpGs: ", format(length(meth_pooled[[cond]]), big.mark=","))
  }
  saveRDS(meth_pooled, meth_cache)
}

# load DMR results
message("\nloading DMR results...")
dmr_results <- list()
for (ct in CONTRASTS) {
  rds_path <- file.path(DMR_DIR, paste0("dmr_", ct$name, ".rds"))
  if (!file.exists(rds_path)) { message("  missing: ", ct$name); next }
  dmr_results[[ct$name]] <- readRDS(rds_path)
  message("  ", ct$name, ": ", length(dmr_results[[ct$name]]), " DMRs")
}

# per-chromosome DMR counts — mirrored bar chart (hypo below, hyper above)
message("\nper-chromosome plots...")
for (ct in CONTRASTS) {
  if (is.null(dmr_results[[ct$name]])) next
  dmrs <- dmr_results[[ct$name]]
  chr_counts <- as.data.frame(table(
    chr=as.character(seqnames(dmrs)), type=dmrs$regionType))
  chr_counts$chr       <- factor(chr_counts$chr, levels=paste0("chr", c(1:22,"X","Y")))
  chr_counts$direction <- ifelse(chr_counts$type=="gain", "Hypo", "Hyper")
  chr_counts$Freq_signed <- ifelse(chr_counts$direction=="Hypo",
                                   -chr_counts$Freq, chr_counts$Freq)
  p <- ggplot(chr_counts, aes(x=chr, y=Freq_signed, fill=direction)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept=0, linewidth=0.4, colour="grey30") +
    scale_fill_manual(values=c(Hypo="#4393C3", Hyper="#D6604D")) +
    scale_y_continuous(labels=function(x) format(abs(x), big.mark=",")) +
    labs(title=paste("DMRs per chromosome --", ct$name),
         x="Chromosome", y="DMRs (hypo below, hyper above)") +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          plot.title=element_text(face="bold"),
          legend.position="top")
  ggsave(file.path(OUT_DIR, paste0(ct$name, "_DMRs_per_chromosome.pdf")),
         p, width=12, height=5)
}

# methylation difference histograms
message("\nmethylation difference plots...")
for (ct in CONTRASTS) {
  if (is.null(dmr_results[[ct$name]])) next
  dmrs <- dmr_results[[ct$name]]
  df <- data.frame(
    methDiff  = dmrs$proportion1 - dmrs$proportion2,
    direction = ifelse(dmrs$regionType=="gain", "Hypo", "Hyper")
  )
  p <- ggplot(df, aes(x=methDiff, fill=direction)) +
    geom_histogram(bins=80, colour="white", linewidth=0.15) +
    geom_vline(xintercept=0, linewidth=0.6, linetype="dashed", colour="grey20") +
    geom_vline(xintercept=c(-0.20, 0.20), linewidth=0.4,
               linetype="dotted", colour="grey50") +
    scale_fill_manual(values=c(Hypo="#4393C3", Hyper="#D6604D")) +
    scale_x_continuous(limits=c(-1,1), breaks=seq(-1,1,0.2),
                       labels=scales::percent_format(accuracy=1)) +
    labs(title=paste("Methylation difference --", ct$name),
         x="Methylation difference (treatment minus reference)",
         y="DMRs") +
    theme_bw(base_size=11) +
    theme(plot.title=element_text(face="bold"), legend.position="top")
  ggsave(file.path(OUT_DIR, paste0(ct$name, "_methylation_difference.pdf")),
         p, width=8, height=5)
}

# locus overlay plots
message("\nlocus plots...")
for (ct in CONTRASTS) {
  if (is.null(dmr_results[[ct$name]])) next
  dmrs   <- dmr_results[[ct$name]]
  meth_a <- meth_pooled[[ct$cond_a]]
  meth_b <- meth_pooled[[ct$cond_b]]
  for (locus in KEY_LOCI) {
    locus_gr   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
    dmrs_locus <- subsetByOverlaps(dmrs, locus_gr)
    out_path   <- file.path(OUT_DIR,
      paste0(ct$name, "_", locus$name, "_locus.pdf"))
    pdf(out_path, width=10, height=5, bg="white")
    par(bg="white", col.axis="black", col.lab="black",
        col.main="black", fg="black")
    tryCatch({
      plotLocalMethylationProfile(
        methylationData1 = meth_a,
        methylationData2 = meth_b,
        region           = locus_gr,
        DMRs             = NULL,
        conditionsNames  = c(ct$cond_a, ct$cond_b),
        windowSize       = 300,
        context          = "CG",
        main             = paste0(locus$name, " -- ", ct$name,
                                  " (", length(dmrs_locus), " DMRs)")
      )
    }, error=function(e) {
      plot.new()
      text(0.5, 0.5, paste0("plot failed: ", conditionMessage(e)),
           cex=0.9, col="red")
    })
    dev.off()
    message("  saved: ", basename(out_path))
  }
}


aso_vpa <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
scr_vpa <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")

# high-confidence filter
hc <- function(gr) gr[mcols(gr)$cytosinesCount >= 6]
aso     <- hc(aso)
aso_vpa <- hc(aso_vpa)
scr_vpa <- hc(scr_vpa)

# build membership matrix for UpSetR
make_membership <- function(sets, names) {
  all_ranges <- reduce(do.call(c, sets))
  mat <- matrix(0L, nrow=length(all_ranges), ncol=length(sets),
                dimnames=list(NULL, names))
  for (i in seq_along(sets)) {
    hits <- findOverlaps(all_ranges, sets[[i]])
    mat[unique(queryHits(hits)), i] <- 1L
  }
  as.data.frame(mat)
}

message("building membership matrix...")
df <- make_membership(
  list(aso, aso_vpa, scr_vpa),
  c("ASO_CTRL", "ASO_VPA", "Scramble_VPA")
)

# main 3-contrast upset plot
pdf(file.path(OUT_DIR, "dmr_upset_plot.pdf"),
    width=12, height=7, onefile=FALSE)
upset(df,
  sets           = c("Scramble_VPA", "ASO_VPA", "ASO_CTRL"),
  keep.order     = TRUE,
  order.by       = "freq",
  sets.bar.color = c("#E67E22", "#C0392B", "#2980B9"),
  main.bar.color = "#2C3E50",
  text.scale     = 1.3,
  mb.ratio       = c(0.55, 0.45),
  mainbar.y.label = "DMR intersections (n)",
  sets.x.label    = "DMRs per contrast",
  point.size     = 3,
  line.size      = 1)
grid::grid.text("DMR overlap across three contrasts (cytosinesCount >= 6)",
  x=0.65, y=0.97,
  gp=grid::gpar(fontsize=11, fontface="bold", col="#1A2A3A"))
dev.off()
message("saved: dmr_upset_plot.pdf")

# summary table
summary_df <- data.frame(
  comparison = c(
    "ASO_CTRL total", "ASO_VPA total", "Scramble_VPA total",
    "ASO_CTRL overlapping ASO_VPA",
    "ASO_CTRL overlapping Scramble_VPA",
    "ASO_VPA overlapping Scramble_VPA",
    "ASO-specific (ASO+ASO_VPA, not Scramble_VPA)"),
  n = c(
    length(aso), length(aso_vpa), length(scr_vpa),
    length(subsetByOverlaps(aso, aso_vpa)),
    length(subsetByOverlaps(aso, scr_vpa)),
    length(subsetByOverlaps(aso_vpa, scr_vpa)),
    length(subsetByOverlaps(
      subsetByOverlaps(aso, aso_vpa), scr_vpa, invert=TRUE)))
)
print(summary_df, row.names=FALSE)
write.table(summary_df, file.path(OUT_DIR, "dmr_overlap_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# supplementary 5-contrast upset plot
message("building 5-contrast matrix...")
aso_ctrl_vpa   <- hc(readRDS("results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds"))
aso_vpa_scrvpa <- hc(readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_VPA.rds"))

df5 <- make_membership(
  list(aso, aso_vpa, scr_vpa, aso_ctrl_vpa, aso_vpa_scrvpa),
  c("ASO_CTRL", "ASO_VPA", "Scramble_VPA", "VPA_on_ASO", "ASO_on_VPA")
)

pdf(file.path(OUT_DIR, "dmr_upset_plot_5contrasts.pdf"),
    width=12, height=8, onefile=FALSE)
upset(df5,
  sets       = c("Scramble_VPA","ASO_VPA","ASO_CTRL","VPA_on_ASO","ASO_on_VPA"),
  keep.order = TRUE,
  order.by   = "freq",
  sets.bar.color = c("#E67E22","#C0392B","#2980B9","#8E44AD","#16A085"),
  main.bar.color = "#2C3E50",
  text.scale     = 1.3,
  mb.ratio       = c(0.6, 0.4),
  mainbar.y.label = "DMR intersections (n)",
  sets.x.label    = "DMRs per contrast")
grid::grid.text("DMR overlap across all 5 contrasts (cytosinesCount >= 6)",
  x=0.65, y=0.985,
  gp=grid::gpar(fontsize=12, fontface="bold", col="#1A2A3A"))
dev.off()
message("saved: dmr_upset_plot_5contrasts.pdf")
message("\ndone. outputs in: ", OUT_DIR)

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

KEEP_CHRS <- paste0("chr", 1:22)
N_BG      <- 1000
set.seed(42)

message("loading ASO-specific DMRs...")
dmrs <- readRDS(DMR_FILE)
dmrs <- dmrs[as.character(seqnames(dmrs)) %in% KEEP_CHRS]
message("  n = ", length(dmrs))

message("loading splice sites...")
bed <- read.table(SPLICE_BED, header=FALSE, sep="\t",
                  col.names=c("chr","start","end"))
bed <- bed[bed$chr %in% KEEP_CHRS, ]
splice_sites <- GRanges(bed$chr, IRanges(bed$end, bed$end))
message("  splice sites: ", length(splice_sites))

# distance from each DMR to nearest splice junction
message("calculating distances...")
dmr_dist      <- distanceToNearest(dmrs, splice_sites)
dmr_distances <- mcols(dmr_dist)$distance
message("  median: ", median(dmr_distances), " bp")
message("  within 300bp: ", sum(dmr_distances <= 300),
        " (", round(100*mean(dmr_distances <= 300), 1), "%)")

# matched background regions
message("generating background (n=", N_BG, ")...")
chr_sizes <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
  chr21=46709983,  chr22=50818468)
chr_counts <- table(as.character(seqnames(dmrs)))
chr_probs  <- chr_counts / sum(chr_counts)
dmr_widths <- width(dmrs)
bg_list  <- list()
attempts <- 0
while (length(bg_list) < N_BG && attempts < N_BG*20) {
  attempts  <- attempts + 1
  chr  <- sample(names(chr_probs), 1, prob=chr_probs)
  w    <- sample(dmr_widths, 1)
  maxs <- chr_sizes[chr] - w
  if (maxs < 1) next
  s <- sample(1:maxs, 1)
  candidate <- GRanges(chr, IRanges(s, s+w-1))
  if (length(findOverlaps(candidate, dmrs)) == 0)
    bg_list[[length(bg_list)+1]] <- candidate
}
bg_regions   <- do.call(c, bg_list)
bg_dist      <- distanceToNearest(bg_regions, splice_sites)
bg_distances <- mcols(bg_dist)$distance
message("  background median: ", median(bg_distances), " bp")

# Wilcoxon test
wtest <- wilcox.test(dmr_distances, bg_distances, alternative="less")
message("Wilcoxon p = ", signif(wtest$p.value, 3))

# save results
results_df <- data.frame(
  type     = c(rep("ASO DMRs",   length(dmr_distances)),
               rep("Background", length(bg_distances))),
  distance = c(dmr_distances, bg_distances)
)
write.csv(results_df, file.path(OUT_DIR, "splice_junction_distances.csv"), row.names=FALSE)

summary_df <- data.frame(
  group            = c("ASO_specific_DMRs", "Background"),
  n                = c(length(dmr_distances), length(bg_distances)),
  median_dist_bp   = c(median(dmr_distances), median(bg_distances)),
  pct_within_300bp = c(round(100*mean(dmr_distances<=300), 1),
                        round(100*mean(bg_distances<=300),  1)),
  wilcox_p         = c(signif(wtest$p.value, 3), NA)
)
print(summary_df, row.names=FALSE)
write.csv(summary_df, file.path(OUT_DIR, "splice_junction_summary.csv"), row.names=FALSE)

# density plot
p <- ggplot(results_df, aes(x=distance, fill=type, colour=type)) +
  geom_density(alpha=0.4, linewidth=0.8) +
  scale_x_log10(labels=scales::comma) +
  scale_fill_manual(values=c("ASO DMRs"="#1B4F8A", "Background"="#cccccc")) +
  scale_colour_manual(values=c("ASO DMRs"="#1B4F8A", "Background"="#888888")) +
  geom_vline(xintercept=300, linetype="dashed", colour="grey40", linewidth=0.5) +
  annotate("text", x=350, y=Inf, label="300bp bin",
           hjust=0, vjust=1.5, size=3, colour="grey40") +
  labs(title="Distance from ASO-specific DMRs to nearest splice junction",
       subtitle=sprintf("n=%d DMRs vs %d background | Wilcoxon p=%s | DMR median=%d bp | BG median=%d bp",
                        length(dmr_distances), N_BG,
                        signif(wtest$p.value,3),
                        median(dmr_distances), median(bg_distances)),
       x="distance to nearest splice junction (bp, log scale)",
       y="density", fill=NULL, colour=NULL) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="top")
ggsave(file.path(OUT_DIR, "splice_junction_distance_density.pdf"), p, width=10, height=6)
message("saved: splice_junction_distance_density.pdf")
message("\ndone. outputs in: ", OUT_DIR)

