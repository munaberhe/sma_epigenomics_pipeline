.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})

setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

# SMN2 locus methylation profile — masked and unmasked alignments
# Uses DMRcaller's plotLocalMethylationProfile for the detailed view
# and computeMethylationProfile for the smoothed lowres overview.
# Replicates are pooled per condition before plotting.

# paths
CHR5_MASKED   <- "results/alignments_smn1_masked/chr5_cx"
BY_CHR_UNMASK <- "results/alignments/bs/by_chr"
OUT_DIR       <- "results/smn2_locus_final"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# parameters
FLANK      <- 2000   # bp either side of the locus for plotting
WIN_SIZE   <- 300    # bin size for detailed profile
WIN_LOWRES <- 500    # bin size for smoothed overview

# SMN1 and SMN2 coordinates on GRCh38
LOCI <- list(
  SMN1 = list(chr="chr5", start=70924941, end=70953015, strand="+"),
  SMN2 = list(chr="chr5", start=70049638, end=70078522, strand="+")
)

# exon coordinates — exon 7 is the ASO target
EXONS <- list(
  SMN1 = data.frame(
    exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start = c(70925030,70938807,70941357,70942326,70942686,
              70944627,70946033,70951913,70952411),
    end   = c(70925158,70938878,70941476,70942526,70942838,
              70944722,70946143,70951966,70952984),
    is_target = c(F,F,F,F,F,F,F,T,F)
  ),
  SMN2 = data.frame(
    exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start = c(70049638,70063415,70065965,70066934,70067294,
              70069235,70070641,70076521,70077019),
    end   = c(70049766,70063486,70066084,70067134,70067446,
              70069330,70070751,70076574,70077592),
    is_target = c(F,F,F,F,F,F,F,T,F)
  )
)

# pairwise comparisons to plot
COMPARISONS <- list(
  list(name="ASO_vs_Scramble_CTRL",
       cond1="ASO_CTRL", cond2="Scramble_CTRL",
       label="ASO effect (CTRL background)"),
  list(name="VPA_vs_CTRL_ASO",
       cond1="ASO_CTRL", cond2="ASO_VPA",
       label="VPA effect (ASO background)"),
  list(name="ASO_vs_Scramble_VPA",
       cond1="ASO_VPA", cond2="Scramble_VPA",
       label="ASO effect (VPA background)")
)

NEEDED <- unique(unlist(lapply(COMPARISONS, function(x) c(x$cond1, x$cond2))))

COND_COLOURS <- c(
  ASO_CTRL      = "#1B4F8A",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#D97706"
)

# build a GRanges with gene and exon features for the plot track
build_gff <- function() {
  rows <- list()
  for (ln in names(LOCI)) {
    l <- LOCI[[ln]]
    rows[[length(rows)+1]] <- data.frame(
      chr=l$chr, start=l$start, end=l$end,
      strand=l$strand, type="gene", name=ln, stringsAsFactors=FALSE)
    ex <- EXONS[[ln]]
    for (i in seq_len(nrow(ex)))
      rows[[length(rows)+1]] <- data.frame(
        chr=l$chr, start=ex$start[i], end=ex$end[i],
        strand=l$strand, type="exon",
        name=sprintf("%s_%s", ln, ex$exon[i]), stringsAsFactors=FALSE)
  }
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          strand=df$strand, type=df$type, name=df$name)
}
GEs <- build_gff()

# read and pool 3 replicates from the SMN1-masked alignment
# filters to CpG context before pooling
read_masked_cpg <- function(condition) {
  message("  masked: ", condition)
  grs <- lapply(1:3, function(r) {
    path <- file.path(CHR5_MASKED,
                      paste0(condition,"_",r,"_chr5.CX_report.txt"))
    d <- read.table(path, header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","countM","countU","context","tri"),
      colClasses=c("character","integer","character","integer",
                   "integer","character","character"))
    d <- d[d$context=="CG", ]
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos,d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

# same but from the original unmasked alignment
read_unmasked_cpg <- function(condition) {
  message("  unmasked: ", condition)
  files <- file.path(BY_CHR_UNMASK,
                     sprintf("%s_%d_chr5.CpG_report.txt.gz", condition, 1:3))
  files <- files[file.exists(files)]
  grs <- lapply(files, function(f) {
    d <- read.table(gzfile(f), header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","countM","countU","context","tri"),
      colClasses=c("character","integer","character","integer",
                   "integer","character","character"))
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos,d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

# detailed locus plot using DMRcaller's plotLocalMethylationProfile
plot_one <- function(pooled, comp, locus_name) {
  locus  <- LOCI[[locus_name]]
  region <- GRanges(seqnames=locus$chr,
                    ranges=IRanges(locus$start-FLANK, locus$end+FLANK))
  plotLocalMethylationProfile(
    methylationData1 = pooled[[comp$cond1]],
    methylationData2 = pooled[[comp$cond2]],
    region           = region,
    DMRs             = NULL,
    conditionsNames  = c(comp$cond1, comp$cond2),
    gff              = GEs,
    windowSize       = WIN_SIZE,
    context          = "CG",
    main             = sprintf("%s: %s vs %s", locus_name, comp$cond1, comp$cond2),
    plotMeanLines    = TRUE,
    plotPoints       = TRUE
  )
  # add exon labels below the gene track
  ex <- EXONS[[locus_name]]
  for (i in seq_len(nrow(ex))) {
    mtext(ex$exon[i], side=1, at=(ex$start[i]+ex$end[i])/2,
          line=-1.5, cex=0.45,
          col=if(ex$is_target[i]) "red" else "black",
          font=if(ex$is_target[i]) 2 else 1)
  }
  usr <- par("usr")
  text(usr[1], usr[3] + (usr[4]-usr[3])*0.05,
       labels=locus_name, cex=0.8, font=2, adj=c(0,0.5))
}

# generate per-comparison PDFs for masked and unmasked
for (alignment in c("masked","unmasked")) {
  message("\n", alignment, " alignment")
  reader <- if(alignment=="masked") read_masked_cpg else read_unmasked_cpg
  pooled <- lapply(NEEDED, reader)
  names(pooled) <- NEEDED

  for (ct in COMPARISONS) {
    fname <- sprintf("SMN_locus_%s_%s.pdf", alignment, ct$name)
    message("  plotting: ", fname)
    pdf(file.path(OUT_DIR, fname), width=11, height=8.5, bg="white")
    par(mfrow=c(2,1), mar=c(5,4,3,1)+0.1, cex=0.9,
        bg="white", col.axis="black", col.lab="black",
        col.main="black", fg="black")
    plot_one(pooled, ct, "SMN1")
    plot_one(pooled, ct, "SMN2")
    dev.off()
  }

  # combined PDF with all comparisons
  fname_all <- sprintf("SMN_locus_%s_all_comparisons.pdf", alignment)
  message("  plotting combined: ", fname_all)
  pdf(file.path(OUT_DIR, fname_all), width=16, height=14, bg="white")
  par(mfrow=c(3,2), mar=c(5,4,3,1)+0.1, cex=0.7,
      bg="white", col.axis="black", col.lab="black",
      col.main="black", fg="black")
  for (ct in COMPARISONS) {
    plot_one(pooled, ct, "SMN1")
    plot_one(pooled, ct, "SMN2")
  }
  dev.off()
}

# weighted mean methylation table for masked data
message("\nweighted means (masked, CpG only)")
masked_pooled <- lapply(NEEDED, read_masked_cpg)
names(masked_pooled) <- NEEDED
rows <- list()
for (cond in NEEDED) {
  for (ln in names(LOCI)) {
    l   <- LOCI[[ln]]
    gr  <- masked_pooled[[cond]]
    sel <- as.character(seqnames(gr))==l$chr &
           start(gr)>=l$start & start(gr)<=l$end & gr$readsN>0
    g   <- gr[sel]
    rows[[length(rows)+1]] <- data.frame(
      condition=cond, locus=ln, n_cpg=length(g),
      weighted_mean=if(length(g)>0) round(sum(g$readsM)/sum(g$readsN),4) else NA)
  }
}
df <- do.call(rbind, rows)
print(df[order(df$locus, df$condition),], row.names=FALSE)
write.table(df, file.path(OUT_DIR, "SMN_weighted_mean_masked.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# smoothed lowres overview using ggplot2 lines (no points)
message("\nlowres smoothed SMN profile")
for (alignment in c("masked","unmasked")) {
  message("  reading ", alignment, " data...")
  reader <- if(alignment=="masked") read_masked_cpg else read_unmasked_cpg
  pooled_lr <- lapply(NEEDED, reader)
  names(pooled_lr) <- NEEDED

  plots <- lapply(names(LOCI), function(ln) {
    locus  <- LOCI[[ln]]
    region <- GRanges(seqnames=locus$chr,
                      ranges=IRanges(locus$start-5000, locus$end+5000))

    # compute binned methylation profile per condition
    prof_df <- do.call(rbind, lapply(NEEDED, function(cond) {
      prof <- computeMethylationProfile(
        methylationData = pooled_lr[[cond]],
        region          = region,
        windowSize      = WIN_LOWRES,
        context         = "CG")
      df           <- as.data.frame(prof)
      df$meth      <- df$sumReadsM / df$sumReadsN
      df$pos       <- (df$start + df$end) / 2
      df$condition <- cond
      df[!is.na(df$meth) & df$sumReadsN >= 3, ]
    }))
    prof_df$condition <- factor(prof_df$condition,
                                levels=c("ASO_CTRL","Scramble_CTRL",
                                         "ASO_VPA","Scramble_VPA"))

    ex7_s <- EXONS[[ln]][EXONS[[ln]]$is_target, "start"]
    ex7_e <- EXONS[[ln]][EXONS[[ln]]$is_target, "end"]

    p <- ggplot(prof_df, aes(x=pos, y=meth, colour=condition)) +
      geom_line(linewidth=0.9, na.rm=TRUE) +
      scale_colour_manual(values=COND_COLOURS) +
      scale_y_continuous(limits=c(0,1), labels=scales::percent_format(1)) +
      scale_x_continuous(labels=function(x) sprintf("%.3f Mb", x/1e6)) +
      theme_classic(base_size=11) +
      theme(legend.position="right",
            panel.grid.major.y=element_line(colour="grey92"),
            plot.title=element_text(face="bold", size=11)) +
      labs(title=ln,
           x=sprintf("Position (chr5, %s alignment)", alignment),
           y="CpG methylation", colour=NULL)

    # mark exon 7 in red for SMN2
    if (ln == "SMN2") {
      p <- p +
        annotate("rect", xmin=ex7_s, xmax=ex7_e, ymin=-Inf, ymax=Inf,
                 fill="red", alpha=0.1) +
        annotate("text", x=(ex7_s+ex7_e)/2, y=0.05,
                 label="E7", colour="red", size=3.5, fontface="bold")
    }
    p
  })

  combined <- (plots[[1]] / plots[[2]]) +
    plot_layout(guides="collect") +
    plot_annotation(
      title    = sprintf("SMN locus methylation — %s alignment", alignment),
      subtitle = sprintf("%dbp bins, pooled replicates, CpG only", WIN_LOWRES),
      theme    = theme(plot.title=element_text(face="bold", size=12),
                       plot.subtitle=element_text(size=9, colour="grey30"))
    ) & theme(legend.position="right")

  fname <- sprintf("SMN_locus_%s_lowres.pdf", alignment)
  ggsave(file.path(OUT_DIR, fname), combined, width=10, height=7)
  message("  saved: ", fname)
}
message("done.")
