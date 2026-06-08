#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

# Low-resolution methylation profiles across all 4 groups.
# Uses the pooled methylation cache from 07b_dmr_plots.R to avoid re-reading raw files.
# Generates: one all-groups plot per region + one pairwise contrast plot per contrast per region.
# Regions: chr1 overview, chrX overview + zoom (ASO DMR hotspot), chr5 + SMN2 zoom.

METH_CACHE <- "results/dmr/meth_pooled_cache.rds"
OUT_DIR    <- "results/lowres_profiles"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

GROUP_COLS <- c(
  ASO_VPA       = "#B2182B",
  ASO_CTRL      = "#1B4F8A",
  Scramble_VPA  = "#F0A500",
  Scramble_CTRL = "#6B7280"
)
GROUP_LTY <- c(
  ASO_VPA       = 1,
  ASO_CTRL      = 1,
  Scramble_VPA  = 2,
  Scramble_CTRL = 2
)
GROUPS <- c("ASO_VPA", "ASO_CTRL", "Scramble_VPA", "Scramble_CTRL")

REGIONS <- list(
  list(region=GRanges("chr1",  IRanges(1, 248956422)),
       window=1000000, label="chr1_1Mb",
       title="CpG methylation chr1 (1 Mb bins)"),
  list(region=GRanges("chrX",  IRanges(1, 155270560)),
       window=500000,  label="chrX_500kb",
       title="CpG methylation chrX (500 kb bins) — ASO hotspot (620 DMRs)"),
  list(region=GRanges("chrX",  IRanges(10000000, 20000000)),
       window=50000,   label="chrX_10_20Mb_50kb",
       title="CpG methylation chrX:10-20 Mb (50 kb bins) — highest ASO density"),
  list(region=GRanges("chr5",  IRanges(1, 181538259)),
       window=500000,  label="chr5_500kb",
       title="CpG methylation chr5 (500 kb bins)"),
  list(region=GRanges("chr5",  IRanges(69000000, 71500000)),
       window=10000,   label="chr5_SMN2_10kb",
       title="CpG methylation chr5:69-71.5 Mb (10 kb bins) — SMN2 locus")
)

CONTRASTS <- list(
  list(g1="ASO_CTRL",     g2="Scramble_CTRL", tag="ASO_effect",
       title="ASO effect (nusinersen vs scramble CTRL)"),
  list(g1="Scramble_VPA", g2="Scramble_CTRL", tag="VPA_effect",
       title="VPA effect (HDAC inhibitor vs CTRL)"),
  list(g1="ASO_VPA",      g2="Scramble_CTRL", tag="combination",
       title="Combination effect (ASO+VPA vs CTRL)"),
  list(g1="ASO_VPA",      g2="ASO_CTRL",      tag="VPA_on_ASO",
       title="VPA effect on ASO background"),
  list(g1="ASO_VPA",      g2="Scramble_VPA",  tag="ASO_on_VPA",
       title="ASO effect on VPA background")
)

message("loading pooled methylation cache...")
if (!file.exists(METH_CACHE)) stop("cache not found: ", METH_CACHE)
meth_pooled <- readRDS(METH_CACHE)
message("  loaded: ", paste(names(meth_pooled), collapse=", "))

for (reg in REGIONS) {
  message("\nregion: ", reg$label)

  # compute profile for each group
  prof_list <- list()
  for (g in GROUPS) {
    message("  computing: ", g)
    prof        <- computeMethylationProfile(meth_pooled[[g]], reg$region,
                                             windowSize=reg$window, context="CG")
    df          <- as.data.frame(prof)
    df$meth     <- df$sumReadsM / df$sumReadsN
    df$pos      <- (df$start + df$end) / 2
    df$group    <- g
    prof_list[[g]] <- df[!is.na(df$meth) & df$sumReadsN >= 3, ]
  }

  x_range <- range(unlist(lapply(prof_list, function(d) d$pos)), na.rm=TRUE)

  # all-groups plot
  out_pdf <- file.path(OUT_DIR, paste0("lowres_allgroups_", reg$label, ".pdf"))
  pdf(out_pdf, width=10, height=4)
  par(mar=c(4,4,3,8)+0.1, bg="white", col.axis="black",
      col.lab="black", col.main="black", fg="black")
  plot(NULL, xlim=x_range, ylim=c(0,1),
       xlab="genomic coordinate", ylab="CpG methylation proportion",
       main=reg$title, cex.main=0.9, font.main=2)
  for (g in GROUPS) {
    d <- prof_list[[g]]
    if (nrow(d) > 0) lines(d$pos, d$meth, col=GROUP_COLS[g], lwd=1.8, lty=GROUP_LTY[g])
  }
  legend("topright", legend=GROUPS, col=GROUP_COLS[GROUPS], lty=GROUP_LTY[GROUPS],
         lwd=1.5, bty="n", cex=0.85, xpd=TRUE, inset=c(-0.15,0))
  dev.off()
  message("  saved: ", basename(out_pdf))

  # pairwise contrast plots
  for (ct in CONTRASTS) {
    out_pdf <- file.path(OUT_DIR,
      paste0("lowres_", ct$tag, "_", reg$label, ".pdf"))
    pdf(out_pdf, width=10, height=4)
    par(mar=c(4,4,3,2)+0.1, bg="white", col.axis="black",
        col.lab="black", col.main="black", fg="black")
    plot(NULL, xlim=x_range, ylim=c(0,1),
         xlab="genomic coordinate", ylab="CpG methylation proportion",
         main=sprintf("%s\n%s (%d kb bins)", ct$title, reg$label, reg$window/1000),
         cex.main=0.9, font.main=2)
    for (g in c(ct$g1, ct$g2)) {
      d <- prof_list[[g]]
      if (nrow(d) > 0) lines(d$pos, d$meth, col=GROUP_COLS[g], lwd=1.8, lty=GROUP_LTY[g])
    }
    legend("topright", legend=c(ct$g1, ct$g2),
           col=GROUP_COLS[c(ct$g1, ct$g2)],
           lty=GROUP_LTY[c(ct$g1, ct$g2)],
           lwd=1.8, bty="n", cex=0.85)
    dev.off()
  }
}


OUT_DIR    <- "results/dmr/plots/annotated"
METH_CACHE <- file.path(DMR_DIR, "meth_pooled_cache.rds")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

WIN_SIZE <- 300

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond_a="Scramble_CTRL", cond_b="ASO_CTRL",
       label="ASO effect (CTRL background)"),
  list(name="ASO_VPA_vs_Scramble_CTRL",
       cond_a="Scramble_CTRL", cond_b="ASO_VPA",
       label="Combination vs baseline"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       cond_a="ASO_CTRL", cond_b="ASO_VPA",
       label="VPA effect on ASO background"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       cond_a="Scramble_VPA", cond_b="ASO_VPA",
       label="ASO effect on VPA background"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       cond_a="Scramble_CTRL", cond_b="Scramble_VPA",
       label="VPA effect (HDAC inhibitor)")
)

LOCI <- list(
  list(name="SLC32A1", chr="chr20", start=38717409, end=38723708,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="SLC32A1 intron (p=5.49e-04, vesicular GABA transporter)",
       annotation="Intron | hypo | +45.5% | neural",
       gene_start=38724486, gene_end=38731000, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="MTA1-DT",   chr="chr14", start=105412884, end=105423739,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="MTA1-DT intron (top hypo hit, p=1.72e-19)",
       annotation="Intron | hypo | +30.8%",
       gene_start=105414884, gene_end=105421739, strand="-",
       exons=data.frame(
         label=c("E1","E2","E3","E4","E5","E6"),
         start=c(105416884,105417581,105417833,105418682,105418939,105419022),
         end=c(105418309,105418309,105418312,105418816,105419080,105419739),
         is_target=rep(FALSE,6))),
  list(name="CHRNB3",    chr="chr8",  start=42702707, end=42709006,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="CHRNB3 exon (p=8.18e-04, nicotinic acetylcholine receptor)",
       annotation="Exon | hypo | +44.2% | neural",
       gene_start=42697376, gene_end=42710436, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="GLRA4",     chr="chrX",  start=103829053, end=103835352,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="GLRA4 chrX (p=5.86e-15, glycine receptor, neural)",
       annotation="Promoter (<=1kb) | hypo | chrX hotspot",
       gene_start=103829053, gene_end=103835352, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="RNF169",    chr="chr11", start=74754377, end=74760676,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="RNF169 intron (p=3.31e-04, top hypo by meth diff +46.3%)",
       annotation="Intron | hypo | +46.3%",
       gene_start=74748849, gene_end=74842414, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="PHACTR3",   chr="chr20", start=59601909,  end=59608208,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="PHACTR3 (p=1.02e-29, top ASO-on-VPA hypo)",
       annotation="Promoter | hypo | ASO_VPA vs Scramble_VPA",
       gene_start=59604540, gene_end=59847711, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="SOX5",      chr="chr12", start=24037319,  end=24043918,
       contrast="ASO_VPA_vs_Scramble_VPA",
       label="SOX5 (p=1.02e-13, neural TF, motor neuron)",
       annotation="Intron | hypo | motor neuron TF",
       gene_start=23533941, gene_end=24090418, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="SEMA3C",    chr="chr7",  start=80810801,  end=80817400,
       contrast="ASO_VPA_vs_Scramble_CTRL",
       label="SEMA3C (p=3.60e-69, axon guidance, neural)",
       annotation="Intron | hypo | axon guidance",
       gene_start=80765501, gene_end=80815949, strand="+",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="GFRA2",     chr="chr8",  start=21810707, end=21817006,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="GFRA2 promoter (p=1.06e-03, GDNF receptor, motor neuron)",
       annotation="Promoter (1-2kb) | hypo | +42.1% | neural",
       gene_start=21692843, gene_end=21812357, strand="-",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0))),
  list(name="SMN2",      chr="chr5",  start=70044638,  end=70083522,
       contrast="ASO_CTRL_vs_Scramble_CTRL",
       label="SMN2 therapeutic target (0 DMRs — null result)",
       annotation="Null result — no DMRs at therapeutic target",
       gene_start=70049638, gene_end=70078522, strand="-",
       exons=data.frame(label=character(0), start=numeric(0),
                        end=numeric(0), is_target=logical(0)))
)

# build GFF for plotLocalMethylationProfile gene track
build_gff <- function(locus) {
  rows <- list()
  rows[[1]] <- data.frame(chr=locus$chr, start=locus$gene_start,
    end=locus$gene_end, strand=locus$strand, type="gene", name=locus$name)
  for (i in seq_len(nrow(locus$exons)))
    rows[[length(rows)+1]] <- data.frame(chr=locus$chr,
      start=locus$exons$start[i], end=locus$exons$end[i],
      strand=locus$strand, type="exon",
      name=paste0(locus$name, "_", locus$exons$label[i]))
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          strand=df$strand, type=df$type, name=df$name)
}

plot_one <- function(meth_a, meth_b, ct, locus, dmrs=NULL) {
  region   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  gff      <- build_gff(locus)
  dmrs_arg <- if (!is.null(dmrs) && length(dmrs) > 0) list("DMRs"=dmrs) else NULL
  COLS <- c(ASO_CTRL="#1B4F8A", ASO_VPA="#B2182B",
            Scramble_VPA="#F0A500", Scramble_CTRL="#6B7280")
  plotLocalMethylationProfile(
    methylationData1 = meth_a,
    methylationData2 = meth_b,
    region           = region,
    DMRs             = dmrs_arg,
    conditionsNames  = c(ct$cond_a, ct$cond_b),
    gff              = gff,
    windowSize       = WIN_SIZE,
    context          = "CG",
    col              = c(COLS[ct$cond_a], COLS[ct$cond_b]),
    main             = sprintf("%s: %s vs %s (%s)\n%s",
                               locus$name, ct$cond_a, ct$cond_b,
                               ct$label, locus$label),
    plotMeanLines = TRUE, plotPoints = TRUE
  )
  # exon labels — staggered to avoid overlap on dense genes
  ex <- locus$exons
  for (i in seq_len(nrow(ex))) {
    mtext(ex$label[i], side=1, at=(ex$start[i]+ex$end[i])/2,
          line=1.2, cex=0.65,
          col=if(ex$is_target[i]) "red" else "black",
          font=if(ex$is_target[i]) 2 else 1)
  }
  usr <- par("usr")
  text(usr[1], usr[3]+(usr[4]-usr[3])*0.05,
       labels=locus$name, cex=0.8, font=2, adj=c(0,0.5))
  # annotation label
  if (!is.null(locus$annotation)) {
    mtext(locus$annotation, side=3, line=0.2, cex=0.7,
          col="grey40", adj=1)
  }
}

message("loading cached pooled data...")
if (!file.exists(METH_CACHE)) stop("cache not found: ", METH_CACHE)
meth_pooled <- readRDS(METH_CACHE)

message("loading DMR results...")
dmr_results <- list()
for (ct in CONTRASTS) {
  rds_path <- file.path(DMR_DIR, paste0("dmr_", ct$name, ".rds"))
  if (!file.exists(rds_path)) { message("  missing: ", ct$name); next }
  dmr_results[[ct$name]] <- readRDS(rds_path)
  message("  ", ct$name, ": ", length(dmr_results[[ct$name]]), " DMRs")
}

# individual PDFs — one per locus using its specific contrast
message("\nindividual locus plots...")
for (locus in LOCI) {
  ct_name <- locus$contrast
  ct <- CONTRASTS[[which(sapply(CONTRASTS, function(x) x$name == ct_name))]]
  if (is.null(ct)) { message("  contrast not found: ", ct_name); next }
  {
    locus_gr   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
    dmrs_locus <- subsetByOverlaps(dmr_results[[ct$name]], locus_gr)
    message("  ", locus$name, " (", ct$name, "): ", length(dmrs_locus), " DMRs in region")
    out_path   <- file.path(OUT_DIR, paste0(locus$name,"_annotated.pdf"))
    cairo_pdf(out_path, width=11, height=6, bg="white")
    par(bg="white", col.axis="black", col.lab="black",
        col.main="black", fg="black", mar=c(5,4,4,2)+0.1)
    tryCatch(
      plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
               ct, locus, dmrs_locus),
      error=function(e) {
        plot.new()
        text(0.5, 0.5, paste0("plot failed: ", conditionMessage(e)),
             col="red", cex=0.9)
      }
    )
    dev.off()
    message("  saved: ", basename(out_path))
  }
}

# combined PDFs — locus-specific contrast only
message("\ncombined per-locus PDFs...")
for (locus in LOCI) {
  ct_name <- locus$contrast
  ct <- CONTRASTS[[which(sapply(CONTRASTS, function(x) x$name == ct_name))]]
  if (is.null(ct)) next
  locus_gr   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  dmrs_locus <- subsetByOverlaps(dmr_results[[ct$name]], locus_gr)
  out_path   <- file.path(OUT_DIR, paste0(locus$name, "_", ct$name, ".pdf"))
  cairo_pdf(out_path, width=11, height=6, bg="white")
  par(bg="white", col.axis="black", col.lab="black",
      col.main="black", fg="black", mar=c(5,4,4,2)+0.1, cex=0.85)
  tryCatch(
    plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
             ct, locus, dmrs_locus),
    error=function(e) { plot.new(); text(0.5,0.5,"failed",col="red") }
  )
  tryCatch(dev.off(), error=function(e) {
    message("  WARNING: dev.off failed for ", basename(out_path), ": ", e$message)
  })
  message("  saved: ", basename(out_path))
}
message("\ndone. outputs in: ", OUT_DIR)
