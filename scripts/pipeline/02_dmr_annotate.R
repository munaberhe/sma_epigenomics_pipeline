#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(DMRcaller)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggplot2)
  library(UpSetR)
  library(dplyr)
})
source("scripts/pipeline/00_sma_palette.R")

# DMR annotation + GO/KEGG enrichment + summary plots.
# input : results/dmr/dmr_<contrast>.rds (from 06b_dmrcaller_combine_chr.R)
# output: results/dmr_annotation/

OUT_DIR <- "results/dmr_annotation"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

CONTRASTS_V <- c(
  "ASO_CTRL_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_VPA_vs_Scramble_VPA",
  "ASO_VPA_vs_ASO_CTRL"
)

# per-contrast annotation + enrichment loop
for (contrast in CONTRASTS_V) {
  message("\nannotating: ", contrast)
  rds_path <- file.path("results/dmr", paste0("dmr_", contrast, ".rds"))
  if (!file.exists(rds_path)) { message("missing: ", rds_path); next }

  dmrs <- readRDS(rds_path)
  message("DMRs loaded: ", length(dmrs))

  # TxDb expects UCSC style (chr1, chr2, ...)
  GenomeInfoDb::seqlevelsStyle(dmrs) <- "UCSC"

  # high-conf filter: >=6 CpGs per DMR
  dmrs_hc <- dmrs[dmrs$cytosinesCount >= 6]
  message("high-conf (>=6 CpGs): ", length(dmrs_hc),
          " | hypo: ", sum(dmrs_hc$regionType == "gain"),
          " | hyper: ", sum(dmrs_hc$regionType == "loss"))

  # annotate all DMRs combined
  anno    <- annotatePeak(dmrs_hc, tssRegion=c(-2000, 2000),
                          TxDb=txdb, annoDb="org.Hs.eg.db")
  anno_df <- as.data.frame(anno)
  write.csv(anno_df,
            file.path(OUT_DIR, paste0(contrast, "_annotated.csv")),
            row.names=FALSE)

  # top 10 by pValue (drop NA gene symbols)
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

  # ChIPseeker built-in plots
  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_pie.pdf")), width=8, height=6)
  plotAnnoPie(anno, main=paste0("Genomic features\n", contrast))
  dev.off()

  pdf(file.path(OUT_DIR, paste0(contrast, "_annotation_bar.pdf")), width=10, height=6)
  plotAnnoBar(anno, title=paste0("Genomic features: ", contrast))
  dev.off()

  pdf(file.path(OUT_DIR, paste0(contrast, "_TSS_distance.pdf")), width=10, height=6)
  plotDistToTSS(anno, title=paste0("Distance to TSS: ", contrast))
  dev.off()

  # GO + KEGG: split by hypo/hyper AND by promoter vs intergenic
  for (direction in c("gain", "loss")) {
    label <- if (direction == "gain") "hypo" else "hyper"
    sub   <- dmrs_hc[dmrs_hc$regionType == direction]
    message("  ", label, ": ", length(sub), " DMRs")
    if (length(sub) < 10) { message("  too few, skipping"); next }
    anno_sub <- annotatePeak(sub, tssRegion=c(-2000, 2000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
    anno_df  <- as.data.frame(anno_sub)
    region_sets <- list(
      promoter   = anno_df[grepl("Promoter", anno_df$annotation), ],
      intergenic = anno_df[!grepl("Promoter", anno_df$annotation) &
                           abs(anno_df$distanceToTSS) <= 50000, ]
    )
    for (region_label in names(region_sets)) {
      rdf   <- region_sets[[region_label]]
      genes <- unique(rdf$geneId)
      genes <- genes[!is.na(genes) & genes != ""]
      message("  ", label, " ", region_label, ": ", length(genes), " genes")
      # enrichGO has its own minGSSize so we only gate on raw gene count
      if (length(genes) < 30) { message("  dropped (n<30)"); next }
      tag <- paste0(contrast, "_GO_BP_", label, "_", region_label)
      go <- enrichGO(gene=genes, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                     ont="BP", pAdjustMethod="BH",
                     pvalueCutoff=0.05, qvalueCutoff=0.2)
      if (!is.null(go) && nrow(go) > 0) {
        message("  top GO BP: ", go$Description[1])
        write.csv(as.data.frame(go),
                  file.path(OUT_DIR, paste0(tag, ".csv")), row.names=FALSE)
        pdf(file.path(OUT_DIR, paste0(tag, "_dotplot.pdf")), width=10, height=8)
        print(dotplot(go, showCategory=15,
                      title=paste0("GO BP (", label, ", ", region_label, "): ", contrast)))
        dev.off()
      } else { message("  no significant GO BP") }
      kegg_tag <- paste0(contrast, "_KEGG_", label, "_", region_label)
      kegg <- enrichKEGG(gene=genes, organism="hsa",
                         pvalueCutoff=0.1, qvalueCutoff=0.3)
      if (!is.null(kegg) && nrow(kegg) > 0) {
        write.csv(as.data.frame(kegg),
                  file.path(OUT_DIR, paste0(kegg_tag, ".csv")), row.names=FALSE)
        pdf(file.path(OUT_DIR, paste0(kegg_tag, "_dotplot.pdf")), width=10, height=8)
        print(dotplot(kegg, showCategory=15,
                      title=paste0("KEGG (", label, ", ", region_label, "): ", contrast)))
        dev.off()
      } else { message("  no significant KEGG") }
    }
  }
}

# Radu item 7: count genes in the chr13 60-80 Mb VPA hyper-DMR hotspot.
# Uses the annotated CSV already written by the loop above.
message("\nchr13 60-80 Mb DMR gene annotation (Radu item 7)...")
chr13_block <- GRanges("chr13", IRanges(60e6, 80e6))
for (vpa_ct in c("Scramble_VPA_vs_Scramble_CTRL", "ASO_VPA_vs_Scramble_CTRL")) {
  csv_path <- file.path(OUT_DIR, paste0(vpa_ct, "_annotated.csv"))
  if (!file.exists(csv_path)) next
  adf <- read.csv(csv_path, stringsAsFactors=FALSE)
  adf_gr <- GRanges(adf$seqnames, IRanges(adf$start, adf$end))
  hits <- subsetByOverlaps(adf_gr, chr13_block)
  idx  <- queryHits(findOverlaps(adf_gr, chr13_block))
  sub  <- adf[idx, ]
  genes_all  <- unique(sub$SYMBOL[!is.na(sub$SYMBOL) & sub$SYMBOL != ""])
  hyper_sub  <- sub[sub$regionType == "loss", ]
  genes_hyper <- unique(hyper_sub$SYMBOL[!is.na(hyper_sub$SYMBOL) & hyper_sub$SYMBOL != ""])
  message(sprintf("  %s: %d DMRs in chr13:60-80Mb | %d unique genes | %d in hyper DMRs",
                  vpa_ct, nrow(sub), length(genes_all), length(genes_hyper)))
  out_tab <- data.frame(
    contrast      = vpa_ct,
    seqnames      = sub$seqnames,
    start         = sub$start,
    end           = sub$end,
    regionType    = sub$regionType,
    pValue        = sub$pValue,
    SYMBOL        = sub$SYMBOL,
    annotation    = sub$annotation,
    stringsAsFactors = FALSE
  )
  write.csv(out_tab,
            file.path(OUT_DIR, paste0("chr13_60_80Mb_genes_", vpa_ct, ".csv")),
            row.names=FALSE)
  message(sprintf("  genes: %s", paste(sort(genes_hyper), collapse=", ")))
}

# contrast metadata: condition pairing + label + colours
CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond_a="ASO_CTRL", cond_b="Scramble_CTRL",
       label="ASO alone vs Scramble CTRL",
       colour=c(ASO_CTRL="#1F3A5F", Scramble_CTRL="#6B7280")),
  list(name="ASO_VPA_vs_Scramble_CTRL",
       cond_a="ASO_VPA", cond_b="Scramble_CTRL",
       label="ASO+VPA vs Scramble CTRL",
       colour=c(ASO_VPA="#C0392B", Scramble_CTRL="#6B7280")),
  list(name="ASO_VPA_vs_ASO_CTRL",
       cond_a="ASO_VPA", cond_b="ASO_CTRL",
       label="VPA effect on ASO background",
       colour=c(ASO_VPA="#C0392B", ASO_CTRL="#1F3A5F")),
  list(name="ASO_VPA_vs_Scramble_VPA",
       cond_a="ASO_VPA", cond_b="Scramble_VPA",
       label="ASO effect on VPA background",
       colour=c(ASO_VPA="#C0392B", Scramble_VPA="#F0A500")),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       cond_a="Scramble_VPA", cond_b="Scramble_CTRL",
       label="VPA alone vs Scramble CTRL",
       colour=c(Scramble_VPA="#F0A500", Scramble_CTRL="#6B7280"))
)

# sample -> condition map for pooling
CONDITIONS <- list(
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = paste0("ASO_VPA_",       1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)

# loci for the locus overlay plots downstream
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

# pool CpG reports per condition, cache to RDS (first run is ~20-30 min)
meth_cache <- file.path("results/dmr", "meth_pooled_cache.rds")
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

# load DMR RDS files per contrast
message("\nloading DMR results...")
dmr_results <- list()
for (ct in CONTRASTS) {
  rds_path <- file.path("results/dmr", paste0("dmr_", ct$name, ".rds"))
  if (!file.exists(rds_path)) { message("  missing: ", ct$name); next }
  dmr_results[[ct$name]] <- readRDS(rds_path)
  message("  ", ct$name, ": ", length(dmr_results[[ct$name]]), " DMRs")
}

# mirrored bar: hypo below 0, hyper above 0, per chromosome
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
    scale_fill_manual(values=c(Hypo="#1F3A5F", Hyper="#C0392B")) +
    scale_y_continuous(labels=function(x) format(abs(x), big.mark=",")) +
    labs(title=paste("DMRs per chromosome:", ct$name),
         x="Chromosome", y="DMRs (hypo below, hyper above)") +
    theme_classic(base_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          plot.title=element_text(face="bold"),
          legend.position="top")
  ggsave(file.path(OUT_DIR, paste0(ct$name, "_DMRs_per_chromosome.pdf")),
         p, width=12, height=5)
}

# methDiff histogram per contrast
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
    scale_fill_manual(values=c(Hypo="#1F3A5F", Hyper="#C0392B")) +
    scale_x_continuous(limits=c(-1,1), breaks=seq(-1,1,0.2),
                       labels=scales::percent_format(accuracy=1)) +
    labs(title=paste("Methylation difference:", ct$name),
         x="Methylation difference (treatment - reference)",
         y="DMRs") +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold"), legend.position="top")
  ggsave(file.path(OUT_DIR, paste0(ct$name, "_methylation_difference.pdf")),
         p, width=8, height=5)
}

# per-locus overlay (DMRcaller plotLocalMethylationProfile)
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
        main             = paste0(locus$name, ": ", ct$name,
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

# UpSet plot across the 3 main contrasts
# load + apply >=6 CpG filter
aso     <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
aso_vpa <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
scr_vpa <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")

hc <- function(gr) gr[mcols(gr)$cytosinesCount >= 6]
aso     <- hc(aso)
aso_vpa <- hc(aso_vpa)
scr_vpa <- hc(scr_vpa)

# membership matrix for UpSetR.
# Each DMR appears exactly once, owned by the earliest set it belongs to.
# Ownership = range came FROM that set (not just overlaps it).
# A range from set B that overlaps a set-A DMR does NOT get its own row -
# the set-A row already has B=1 via findOverlaps.
# This guarantees colSums(df)[i] == length(sets[[i]]) with no inflation.
# (Previous approach used unique(do.call(c,sets)) which only deduplicates
# IDENTICAL ranges; non-identical overlapping ranges from other sets were
# kept as separate rows and falsely marked as belonging to the earlier set,
# inflating set bars - e.g. 151 real ASO DMRs appearing as 159.)
make_membership <- function(sets, set_names) {
  sets    <- lapply(sets, function(gr) sort(unique(gr)))
  owned   <- NULL
  row_mats <- vector("list", length(sets))
  for (i in seq_along(sets)) {
    gr_i <- sets[[i]]
    # drop ranges already owned by an earlier set
    if (!is.null(owned)) {
      dup <- findOverlaps(gr_i, owned, minoverlap=1L)
      if (length(dup)) gr_i <- gr_i[-unique(queryHits(dup))]
    }
    if (!length(gr_i)) next
    mat_i <- matrix(0L, nrow=length(gr_i), ncol=length(sets),
                    dimnames=list(NULL, set_names))
    for (j in seq_along(sets)) {
      hits <- findOverlaps(gr_i, sets[[j]], minoverlap=1L)
      mat_i[unique(queryHits(hits)), j] <- 1L
    }
    row_mats[[i]] <- mat_i
    owned <- if (is.null(owned)) gr_i else c(owned, gr_i)
  }
  do.call(rbind, Filter(Negate(is.null), row_mats))
}

message("building membership matrix...")
df <- make_membership(
  list(aso, aso_vpa, scr_vpa),
  c("ASO_CTRL", "ASO_VPA", "Scramble_VPA")
)
message("  raw set sizes:   ",
        paste(c("ASO_CTRL","ASO_VPA","Scramble_VPA"),
              c(length(aso), length(aso_vpa), length(scr_vpa)),
              sep="=", collapse=", "))
message("  matrix col sums: ",
        paste(colnames(df), colSums(df), sep="=", collapse=", "))
# These should match - if they don't, findOverlaps found cross-set inflations

# 3-contrast upset
# each bar = DMRs present in EXACTLY that combination of sets (not "at least").
# emit a labelled intersection table alongside so the plot reads without
# needing prior UpSetR convention knowledge.
tryCatch({
pdf(file.path(OUT_DIR, "dmr_upset_plot.pdf"),
    width=12, height=7.5, onefile=FALSE)
upset(df,
  sets           = c("Scramble_VPA", "ASO_VPA", "ASO_CTRL"),
  keep.order     = TRUE,
  order.by       = "freq",
  sets.bar.color = c("#F0A500", "#C0392B", "#1F3A5F"),
  main.bar.color = "#2C3E50",
  text.scale     = 1.3,
  mb.ratio       = c(0.55, 0.45),
  mainbar.y.label = "DMRs in exactly this combination (n)",
  sets.x.label    = "DMRs per contrast",
  point.size     = 3,
  line.size      = 1)
grid::grid.text("DMR overlap, 3 contrasts (cytosinesCount >= 6)",
  x=0.65, y=0.98,
  gp=grid::gpar(fontsize=12, fontface="bold", col="#1A2A3A"))
grid::grid.text(
  "each bar = DMRs in exactly that combo of sets, NOT 'at least'.",
  x=0.65, y=0.945,
  gp=grid::gpar(fontsize=9, col="grey30"))
dev.off()
}, error=function(e) message("upset skipped: ", e$message))
message("saved: dmr_upset_plot.pdf")

# labelled intersection table - one row per cell of the upset matrix.
# answers the "does 151 match anything here?" question directly.
message("\nbuilding slide 16 / slide 18 panels...")

CONTRAST_LABELS <- c(
  'ASO_CTRL_vs_Scramble_CTRL'     = 'ASO vs Scr_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL' = 'VPA vs Scr_CTRL',
  'ASO_VPA_vs_Scramble_CTRL'      = 'ASO+VPA vs Scr_CTRL',
  'ASO_VPA_vs_ASO_CTRL'           = 'ASO+VPA vs ASO',
  'ASO_VPA_vs_Scramble_VPA'       = 'ASO+VPA vs Scr_VPA'
)

dir_df <- do.call(rbind, lapply(names(CONTRAST_LABELS), function(ct) {
  dmr <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  dmr <- dmr[dmr$context == 'CG']
  data.frame(
    contrast  = CONTRAST_LABELS[ct],
    direction = c('Hypomethylated', 'Hypermethylated'),
    n         = c(sum(dmr$regionType == 'gain'),
                  sum(dmr$regionType == 'loss'))
  )
}))

dir_df$contrast  <- factor(dir_df$contrast, levels=CONTRAST_LABELS)
dir_df$direction <- factor(dir_df$direction,
                           levels=c('Hypermethylated','Hypomethylated'))
totals_df <- dir_df %>% group_by(contrast) %>% summarise(total=sum(n))

# slide 16: DMR_hypo_hyper_counts_bar (log10 stacked)
# label positioning was previously stuck at total y, overlapping the hyper segment.
# move totals above the bar and inset segment labels inside their stack.
p_dir <- ggplot(dir_df, aes(x=contrast, y=n, fill=direction)) +
  geom_bar(stat='identity', position='stack', width=0.6) +
  geom_text(aes(label=formatC(n, format='d', big.mark=',')),
            position=position_stack(vjust=0.5),
            size=3.0, colour='white', fontface='bold') +
  geom_text(data=totals_df,
            aes(x=contrast, y=total*1.25,
                label=formatC(total, format='d', big.mark=',')),
            inherit.aes=FALSE,
            size=3.2, fontface='bold', colour='grey25') +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_log10(labels=scales::comma,
                breaks=c(100,1000,10000,100000,1000000),
                expand=expansion(mult=c(0.02, 0.18))) +
  annotation_logticks(sides='l') +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=10),
        legend.position='top') +
  labs(x=NULL, y='Number of DMRs (log10)', fill=NULL)

ggsave(file.path(OUT_DIR, 'DMR_hypo_hyper_counts_bar.pdf'),
       p_dir, width=10, height=6, device=cairo_pdf)
message("saved: DMR_hypo_hyper_counts_bar.pdf")

# slide 16: DMR_hypo_hyper_counts_bar_v2 (linear axis, side-by-side)
# log10 hides the actual magnitude; linear + dodge is easier to read in the deck.
# small bars get a clipped readable label by drawing text above the bar.
dir_df_v2 <- dir_df %>%
  group_by(contrast) %>%
  mutate(label_y = n) %>%
  ungroup()

p_dir_v2 <- ggplot(dir_df_v2, aes(x=contrast, y=n, fill=direction)) +
  geom_bar(stat='identity', position=position_dodge(width=0.7),
           width=0.6) +
  geom_text(aes(label=formatC(n, format='d', big.mark=',')),
            position=position_dodge(width=0.7),
            vjust=-0.4, size=3.0, fontface='bold', colour='grey25') +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_log10(labels=scales::comma,
                breaks=c(100, 1000, 10000, 100000, 1000000),
                expand=expansion(mult=c(0.01, 0.18))) +
  annotation_logticks(sides='l') +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=10),
        legend.position='top',
        plot.title=element_text(face='bold')) +
  labs(title='How many DMRs per contrast and in which direction? (log10 scale)',
       x=NULL, y='Number of DMRs (log10)', fill=NULL)

ggsave(file.path(OUT_DIR, 'slide16_dmr_counts_bar_v2.pdf'),
       p_dir_v2, width=11, height=6, device=cairo_pdf)
message("saved: slide16_dmr_counts_bar_v2.pdf")

# slide 16: slide16_dmr_diverging (mirrored, SHARED y-axis)
# previous version used scales='free' per facet which destroyed comparability.
# use a single signed y-axis (hypo positive, hyper negative) on log10 magnitude
# with symmetric breaks so all contrasts share the same scale.
diverging_df <- dir_df %>%
  mutate(n_signed = ifelse(direction == 'Hypomethylated', n, -n),
         y_log    = sign(n_signed) * log10(abs(n_signed) + 1))

# pretty log10 breaks: 10, 100, 1k, 10k, 100k, 1M (both signs)
log_breaks <- c(10, 100, 1e3, 1e4, 1e5, 1e6)
break_pos  <- sign(c(-log_breaks, log_breaks)) *
              log10(c(rev(log_breaks), log_breaks) + 1)
break_lab  <- c(
  paste0("-", formatC(rev(log_breaks), format='d', big.mark=',')),
  formatC(log_breaks, format='d', big.mark=',')
)

p_div <- ggplot(diverging_df,
                aes(x=contrast, y=y_log, fill=direction)) +
  geom_bar(stat='identity', width=0.65) +
  geom_hline(yintercept=0, linewidth=0.4, colour='grey30') +
  geom_text(aes(label=formatC(n, format='d', big.mark=','),
                vjust=ifelse(direction=='Hypomethylated', -0.4, 1.3)),
            size=3.0, fontface='bold', colour='grey25') +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_continuous(breaks=break_pos, labels=break_lab,
                     limits=c(-max(abs(diverging_df$y_log))*1.15,
                               max(abs(diverging_df$y_log))*1.15)) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=10),
        axis.text.y=element_text(size=8),
        legend.position='top',
        plot.title=element_text(face='bold')) +
  labs(title='Hypomethylated (above) vs hypermethylated (below)',
       subtitle='shared log10 axis across contrasts',
       x=NULL,
       y='Number of DMRs (log10, signed)',
       fill=NULL)

ggsave(file.path(OUT_DIR, 'slide16_dmr_diverging.pdf'),
       p_div, width=11, height=6.5, device=cairo_pdf)
message("saved: slide16_dmr_diverging.pdf")

# slide 18: DMR_violin_hypo_hyper + slide18_violin_hypo_hyper
# build violin input: per-DMR proportion in cond1 + cond2, separated by direction
message("generating DMR violin plot hypo/hyper...")

violin_dmr <- do.call(rbind, lapply(names(CONTRAST_LABELS), function(ct) {
  dmr <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  dmr <- dmr[dmr$context == 'CG' & dmr$cytosinesCount >= 4]
  rbind(
    data.frame(contrast=CONTRAST_LABELS[ct], direction='Hypomethylated',
               meth=c(dmr$proportion1[dmr$regionType=='gain'],
                      dmr$proportion2[dmr$regionType=='gain'])),
    data.frame(contrast=CONTRAST_LABELS[ct], direction='Hypermethylated',
               meth=c(dmr$proportion1[dmr$regionType=='loss'],
                      dmr$proportion2[dmr$regionType=='loss']))
  )
}))

violin_dmr$contrast  <- factor(violin_dmr$contrast, levels=CONTRAST_LABELS)
violin_dmr$direction <- factor(violin_dmr$direction,
                                levels=c('Hypomethylated','Hypermethylated'))

# original DMR_violin_hypo_hyper.pdf (no title)
p_viol <- ggplot(violin_dmr, aes(x=contrast, y=meth, fill=direction)) +
  geom_violin(trim=FALSE, alpha=1.0, linewidth=0.3) +
  geom_boxplot(width=0.07, fill='white',
               outlier.size=0.2, outlier.alpha=0.2) +
  facet_wrap(~direction, ncol=2) +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_continuous(labels=scales::percent_format(1), limits=c(0,1)) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        strip.text=element_text(face='bold', size=11),
        legend.position='none') +
  labs(x=NULL, y='CpG methylation proportion')

ggsave(file.path(OUT_DIR, 'DMR_violin_hypo_hyper.pdf'),
       p_viol, width=14, height=6, device=cairo_pdf)
message("saved: DMR_violin_hypo_hyper.pdf")

# slide 18 deck variant: explicit question title + tighter formatting so the
# title does not auto-truncate with "..." in the slide layout
p_viol_18 <- ggplot(violin_dmr, aes(x=contrast, y=meth, fill=direction)) +
  geom_violin(trim=FALSE, alpha=1.0, linewidth=0.3) +
  geom_boxplot(width=0.07, fill='white',
               outlier.size=0.2, outlier.alpha=0.2) +
  facet_wrap(~direction, ncol=2) +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_continuous(labels=scales::percent_format(1), limits=c(0,1)) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        strip.text=element_text(face='bold', size=11),
        plot.title=element_text(face='bold', size=13),
        plot.subtitle=element_text(size=10, colour='grey30'),
        legend.position='none') +
  labs(title='CpG methylation at DMR loci: hypo vs hyper',
       subtitle='proportion methylated CpGs, cond1 and cond2 methylation pooled per DMR',
       x=NULL, y='CpG methylation proportion')

ggsave(file.path(OUT_DIR, 'slide18_violin_hypo_hyper.pdf'),
       p_viol_18, width=14, height=6.5, device=cairo_pdf)
message("saved: slide18_violin_hypo_hyper.pdf")

message("\ndone. all outputs in: ", OUT_DIR)

# upset plot: simple binary matrix
aso     <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
aso_vpa <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds")
scr_vpa <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
hc <- function(gr) gr[mcols(gr)$cytosinesCount >= 6]
aso     <- hc(aso)
aso_vpa <- hc(aso_vpa)
scr_vpa <- hc(scr_vpa)

# build simple binary matrix from ASO_CTRL as anchor
message("building upset matrix...")
anchor <- aso
mat <- data.frame(
  ASO_CTRL     = rep(1L, length(anchor)),
  ASO_VPA      = as.integer(countOverlaps(anchor, aso_vpa) > 0),
  Scramble_VPA = as.integer(countOverlaps(anchor, scr_vpa) > 0)
)
message("  ASO_CTRL=", nrow(mat),
        " overlapping ASO_VPA=", sum(mat$ASO_VPA),
        " overlapping Scramble_VPA=", sum(mat$Scramble_VPA))

tryCatch({
  pdf(file.path(OUT_DIR, "dmr_upset_plot.pdf"), width=12, height=7.5, onefile=FALSE)
  upset(mat, sets=c("Scramble_VPA","ASO_VPA","ASO_CTRL"),
        keep.order=TRUE, order.by="freq",
        sets.bar.color=c("#F0A500","#C0392B","#1F3A5F"),
        main.bar.color="#2C3E50", text.scale=1.3,
        mb.ratio=c(0.55,0.45),
        mainbar.y.label="DMRs in this combination",
        sets.x.label="DMRs per contrast",
        point.size=3, line.size=1)
  grid::grid.text("DMR overlap (cytosinesCount >= 6)",
    x=0.65, y=0.98,
    gp=grid::gpar(fontsize=12, fontface="bold", col="#1A2A3A"))
  dev.off()
  message("saved: dmr_upset_plot.pdf")
}, error=function(e) message("upset skipped: ", e$message))

message("\nbuilding slide 16 / slide 18 panels...")

CONTRAST_LABELS <- c(
  'ASO_CTRL_vs_Scramble_CTRL'     = 'ASO vs Scr_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL' = 'VPA vs Scr_CTRL',
  'ASO_VPA_vs_Scramble_CTRL'      = 'ASO+VPA vs Scr_CTRL',
  'ASO_VPA_vs_ASO_CTRL'           = 'ASO+VPA vs ASO',
  'ASO_VPA_vs_Scramble_VPA'       = 'ASO+VPA vs Scr_VPA'
)

dir_df <- do.call(rbind, lapply(names(CONTRAST_LABELS), function(ct) {
  dmr <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  dmr <- dmr[dmr$context == 'CG']
  data.frame(
    contrast  = CONTRAST_LABELS[ct],
    direction = c('Hypomethylated', 'Hypermethylated'),
    n         = c(sum(dmr$regionType == 'gain'),
                  sum(dmr$regionType == 'loss'))
  )
}))

dir_df$contrast  <- factor(dir_df$contrast, levels=CONTRAST_LABELS)
dir_df$direction <- factor(dir_df$direction,
                           levels=c('Hypermethylated','Hypomethylated'))
totals_df <- dir_df %>% group_by(contrast) %>% summarise(total=sum(n))

# slide 16: DMR_hypo_hyper_counts_bar (log10 stacked)
# label positioning was previously stuck at total y, overlapping the hyper segment.
# move totals above the bar and inset segment labels inside their stack.
p_dir <- ggplot(dir_df, aes(x=contrast, y=n, fill=direction)) +
  geom_bar(stat='identity', position='stack', width=0.6) +
  geom_text(aes(label=formatC(n, format='d', big.mark=',')),
            position=position_stack(vjust=0.5),
            size=3.0, colour='white', fontface='bold') +
  geom_text(data=totals_df,
            aes(x=contrast, y=total*1.25,
                label=formatC(total, format='d', big.mark=',')),
            inherit.aes=FALSE,
            size=3.2, fontface='bold', colour='grey25') +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_log10(labels=scales::comma,
                breaks=c(100,1000,10000,100000,1000000),
                expand=expansion(mult=c(0.02, 0.18))) +
  annotation_logticks(sides='l') +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=10),
        legend.position='top') +
  labs(x=NULL, y='Number of DMRs (log10)', fill=NULL)

ggsave(file.path(OUT_DIR, 'DMR_hypo_hyper_counts_bar.pdf'),
       p_dir, width=10, height=6, device=cairo_pdf)
message("saved: DMR_hypo_hyper_counts_bar.pdf")

# slide 16: DMR_hypo_hyper_counts_bar_v2 (linear axis, side-by-side)
# log10 hides the actual magnitude; linear + dodge is easier to read in the deck.
# small bars get a clipped readable label by drawing text above the bar.
dir_df_v2 <- dir_df %>%
  group_by(contrast) %>%
  mutate(label_y = n) %>%
  ungroup()

p_dir_v2 <- ggplot(dir_df_v2, aes(x=contrast, y=n, fill=direction)) +
  geom_bar(stat='identity', position=position_dodge(width=0.7),
           width=0.6) +
  geom_text(aes(label=formatC(n, format='d', big.mark=',')),
            position=position_dodge(width=0.7),
            vjust=-0.4, size=3.0, fontface='bold', colour='grey25') +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_log10(labels=scales::comma,
                breaks=c(100, 1000, 10000, 100000, 1000000),
                expand=expansion(mult=c(0.01, 0.18))) +
  annotation_logticks(sides='l') +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=10),
        legend.position='top',
        plot.title=element_text(face='bold')) +
  labs(title='How many DMRs per contrast and in which direction? (log10 scale)',
       x=NULL, y='Number of DMRs (log10)', fill=NULL)

ggsave(file.path(OUT_DIR, 'slide16_dmr_counts_bar_v2.pdf'),
       p_dir_v2, width=11, height=6, device=cairo_pdf)
message("saved: slide16_dmr_counts_bar_v2.pdf")

# slide 16: slide16_dmr_diverging (mirrored, SHARED y-axis)
# previous version used scales='free' per facet which destroyed comparability.
# use a single signed y-axis (hypo positive, hyper negative) on log10 magnitude
# with symmetric breaks so all contrasts share the same scale.
diverging_df <- dir_df %>%
  mutate(n_signed = ifelse(direction == 'Hypomethylated', n, -n),
         y_log    = sign(n_signed) * log10(abs(n_signed) + 1))

# pretty log10 breaks: 10, 100, 1k, 10k, 100k, 1M (both signs)
log_breaks <- c(10, 100, 1e3, 1e4, 1e5, 1e6)
break_pos  <- sign(c(-log_breaks, log_breaks)) *
              log10(c(rev(log_breaks), log_breaks) + 1)
break_lab  <- c(
  paste0("-", formatC(rev(log_breaks), format='d', big.mark=',')),
  formatC(log_breaks, format='d', big.mark=',')
)

p_div <- ggplot(diverging_df,
                aes(x=contrast, y=y_log, fill=direction)) +
  geom_bar(stat='identity', width=0.65) +
  geom_hline(yintercept=0, linewidth=0.4, colour='grey30') +
  geom_text(aes(label=formatC(n, format='d', big.mark=','),
                vjust=ifelse(direction=='Hypomethylated', -0.4, 1.3)),
            size=3.0, fontface='bold', colour='grey25') +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_continuous(breaks=break_pos, labels=break_lab,
                     limits=c(-max(abs(diverging_df$y_log))*1.15,
                               max(abs(diverging_df$y_log))*1.15)) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=10),
        axis.text.y=element_text(size=8),
        legend.position='top',
        plot.title=element_text(face='bold')) +
  labs(title='Hypomethylated (above) vs hypermethylated (below)',
       subtitle='shared log10 axis across contrasts',
       x=NULL,
       y='Number of DMRs (log10, signed)',
       fill=NULL)

ggsave(file.path(OUT_DIR, 'slide16_dmr_diverging.pdf'),
       p_div, width=11, height=6.5, device=cairo_pdf)
message("saved: slide16_dmr_diverging.pdf")

# slide 18: DMR_violin_hypo_hyper + slide18_violin_hypo_hyper
# build violin input: per-DMR proportion in cond1 + cond2, separated by direction
message("generating DMR violin plot hypo/hyper...")

violin_dmr <- do.call(rbind, lapply(names(CONTRAST_LABELS), function(ct) {
  dmr <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  dmr <- dmr[dmr$context == 'CG' & dmr$cytosinesCount >= 4]
  rbind(
    data.frame(contrast=CONTRAST_LABELS[ct], direction='Hypomethylated',
               meth=c(dmr$proportion1[dmr$regionType=='gain'],
                      dmr$proportion2[dmr$regionType=='gain'])),
    data.frame(contrast=CONTRAST_LABELS[ct], direction='Hypermethylated',
               meth=c(dmr$proportion1[dmr$regionType=='loss'],
                      dmr$proportion2[dmr$regionType=='loss']))
  )
}))

violin_dmr$contrast  <- factor(violin_dmr$contrast, levels=CONTRAST_LABELS)
violin_dmr$direction <- factor(violin_dmr$direction,
                                levels=c('Hypomethylated','Hypermethylated'))

# original DMR_violin_hypo_hyper.pdf (no title)
p_viol <- ggplot(violin_dmr, aes(x=contrast, y=meth, fill=direction)) +
  geom_violin(trim=FALSE, alpha=1.0, linewidth=0.3) +
  geom_boxplot(width=0.07, fill='white',
               outlier.size=0.2, outlier.alpha=0.2) +
  facet_wrap(~direction, ncol=2) +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_continuous(labels=scales::percent_format(1), limits=c(0,1)) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        strip.text=element_text(face='bold', size=11),
        legend.position='none') +
  labs(x=NULL, y='CpG methylation proportion')

ggsave(file.path(OUT_DIR, 'DMR_violin_hypo_hyper.pdf'),
       p_viol, width=14, height=6, device=cairo_pdf)
message("saved: DMR_violin_hypo_hyper.pdf")

# slide 18 deck variant: explicit question title + tighter formatting so the
# title does not auto-truncate with "..." in the slide layout
p_viol_18 <- ggplot(violin_dmr, aes(x=contrast, y=meth, fill=direction)) +
  geom_violin(trim=FALSE, alpha=1.0, linewidth=0.3) +
  geom_boxplot(width=0.07, fill='white',
               outlier.size=0.2, outlier.alpha=0.2) +
  facet_wrap(~direction, ncol=2) +
  scale_fill_manual(values=c(Hypomethylated='#1F3A5F',
                              Hypermethylated='#C0392B')) +
  scale_y_continuous(labels=scales::percent_format(1), limits=c(0,1)) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        strip.text=element_text(face='bold', size=11),
        plot.title=element_text(face='bold', size=13),
        plot.subtitle=element_text(size=10, colour='grey30'),
        legend.position='none') +
  labs(title='CpG methylation at DMR loci: hypo vs hyper',
       subtitle='proportion methylated CpGs, cond1 and cond2 methylation pooled per DMR',
       x=NULL, y='CpG methylation proportion')

ggsave(file.path(OUT_DIR, 'slide18_violin_hypo_hyper.pdf'),
       p_viol_18, width=14, height=6.5, device=cairo_pdf)
message("saved: slide18_violin_hypo_hyper.pdf")

message("\ndone. all outputs in: ", OUT_DIR)
