#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
})

# DMR diagnostic and locus plots.
# Generates per-chromosome bar charts, methylation difference histograms,
# and locus overlay plots using plotLocalMethylationProfile.

DMR_DIR    <- "results/dmr"
BY_CHR_DIR <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr/plots"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

CHROMS <- paste0("chr", c(1:22, "X", "Y"))

# DMRcaller direction note: "gain" = treatment hypomethylated, "loss" = hypermethylated
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
message("\ndone. outputs in: ", OUT_DIR)
