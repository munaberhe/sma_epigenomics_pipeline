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
  list(name="VPA_vs_Scramble_CTRL",
       cond1="Scramble_VPA", cond2="Scramble_CTRL",
       label="VPA effect (CTRL background)"),
  list(name="ASO_VPA_vs_Scramble_CTRL",
       cond1="ASO_VPA", cond2="Scramble_CTRL",
       label="Combined vs CTRL"),
  list(name="VPA_vs_CTRL_ASO",
       cond1="ASO_CTRL", cond2="ASO_VPA",
       label="VPA effect (ASO background)"),
  list(name="ASO_vs_Scramble_VPA",
       cond1="ASO_VPA", cond2="Scramble_VPA",
       label="ASO effect (VPA background)")
)

NEEDED <- unique(unlist(lapply(COMPARISONS, function(x) c(x$cond1, x$cond2))))

SENSITIVE_DMRS <- list(
  ASO_vs_Scramble_VPA = GRanges(
    seqnames   = "chr5",
    ranges     = IRanges(start=c(70074938, 70088438), end=c(70074987, 70088487)),
    regionType = c("loss", "gain"),
    pValue     = c(0.037, 0.033)
  ),
  VPA_vs_Scramble_CTRL = GRanges(
    seqnames   = "chr5",
    ranges     = IRanges(start=c(70046438, 70074838, 70079338),
                         end  =c(70046487, 70074887, 70079387)),
    regionType = c("gain", "gain", "gain"),
    pValue     = c(0.011, 0.002, 0.011)
  )
)

COND_COLOURS <- c(
  ASO_CTRL      = "#1B4F8A",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#B2182B",
  Scramble_VPA  = "#F0A500"
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
plot_one <- function(pooled, comp, locus_name, dmrs=NULL) {
  locus  <- LOCI[[locus_name]]
  region <- GRanges(seqnames=locus$chr,
                    ranges=IRanges(locus$start-FLANK, locus$end+FLANK))
  plotLocalMethylationProfile(
    methylationData1 = pooled[[comp$cond1]],
    methylationData2 = pooled[[comp$cond2]],
    region           = region,
    DMRs             = if(!is.null(dmrs)) list("DMRs"=dmrs) else NULL,
    conditionsNames  = c(comp$cond1, comp$cond2),
    gff              = GEs,
    windowSize       = WIN_SIZE,
    context          = "CG",
    col              = c(COND_COLOURS[comp$cond1], COND_COLOURS[comp$cond2],
                           "#FFFFFF", "#FFFFFF", "#E41A1C"),
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
    h <- if(alignment=="masked") 6 else 8.5
    nr <- if(alignment=="masked") 1 else 2
    cairo_pdf(file.path(OUT_DIR, fname), width=11, height=h)
    par(mfrow=c(nr,1), mar=c(5,4,3,1)+0.1, cex=0.9,
        bg="white", col.axis="black", col.lab="black",
        col.main="black", fg="black")
    if(alignment=="unmasked") plot_one(pooled, ct, "SMN1")
    ct_key <- gsub("ASO_VPA_vs_Scramble_VPA", "ASO_vs_Scramble_VPA",
                   gsub("Scramble_VPA_vs_Scramble_CTRL", "VPA_vs_Scramble_CTRL", ct$name))
    dmrs_arg <- if(!is.null(SENSITIVE_DMRS[[ct_key]])) SENSITIVE_DMRS[[ct_key]] else NULL
    plot_one(pooled, ct, "SMN2", dmrs=dmrs_arg)
    dev.off()
  }

  # combined PDF with all comparisons
  fname_all <- sprintf("SMN_locus_%s_all_comparisons.pdf", alignment)
  message("  plotting combined: ", fname_all)
  n_loci   <- if(alignment=="masked") 1 else 2
  n_panels <- length(COMPARISONS) * n_loci
  n_cols <- 2
  n_rows <- ceiling(n_panels / n_cols)
  cairo_pdf(file.path(OUT_DIR, fname_all), width=11, height=n_rows*3.5)
  par(mfrow=c(n_rows, n_cols), mar=c(5,4,3,1)+0.1, cex=0.7,
      bg="white", col.axis="black", col.lab="black",
      col.main="black", fg="black")
  for (ct in COMPARISONS) {
    if(alignment=="unmasked") plot_one(pooled, ct, "SMN1")
    ct_key <- gsub("ASO_VPA_vs_Scramble_VPA", "ASO_vs_Scramble_VPA",
                   gsub("Scramble_VPA_vs_Scramble_CTRL", "VPA_vs_Scramble_CTRL", ct$name))
    dmrs_arg <- if(!is.null(SENSITIVE_DMRS[[ct_key]])) SENSITIVE_DMRS[[ct_key]] else NULL
    plot_one(pooled, ct, "SMN2", dmrs=dmrs_arg)
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

df_both <- rbind(df_masked, df_unmasked)
df_both$condition <- factor(df_both$condition,
  levels=c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA"))

make_panel <- function(df, title) {
  ex7_s <- EXONS[["SMN2"]][EXONS[["SMN2"]]$is_target, "start"]
  ex7_e <- EXONS[["SMN2"]][EXONS[["SMN2"]]$is_target, "end"]
  ggplot(df, aes(x=pos, y=meth, colour=condition)) +
    geom_line(linewidth=0.9, na.rm=TRUE) +
    annotate("rect", xmin=ex7_s, xmax=ex7_e, ymin=-Inf, ymax=Inf,
             fill="red", alpha=0.1) +
    annotate("text", x=(ex7_s+ex7_e)/2, y=0.05,
             label="E7", colour="red", size=3.5, fontface="bold") +
    scale_colour_manual(values=COND_COLOURS) +
    scale_y_continuous(limits=c(0,1), labels=scales::percent_format(1)) +
    scale_x_continuous(labels=function(x) sprintf("%.3f Mb", x/1e6)) +
    theme_classic(base_size=11) +
    theme(legend.position="right",
          panel.grid.major.y=element_line(colour="grey92"),
          plot.title=element_text(face="bold", size=11)) +
    labs(title=title, x="Position (chr5)", y="CpG methylation", colour=NULL)
}

p_masked   <- make_panel(df_masked,   "SMN2 — SMN1-masked alignment")
p_unmasked <- make_panel(df_unmasked, "SMN2 — unmasked alignment")

combined <- (p_unmasked / p_masked) +
  plot_layout(guides="collect") &
  theme(legend.position="right")

ggsave(file.path(OUT_DIR, "SMN2_masked_vs_unmasked_comparison.pdf"),
       combined, width=12, height=8, device=cairo_pdf)
message("  saved: SMN2_masked_vs_unmasked_comparison.pdf")
