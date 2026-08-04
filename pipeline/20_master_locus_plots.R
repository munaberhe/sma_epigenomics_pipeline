#!/usr/bin/env Rscript
# ============================================================================
# 20_master_locus_plots.R
# Single script for ALL locus plots - SMN2, SMN2 sensitive, survivors,
# SEMA3C/SIGMAR1, and new pairwise candidates.
# One plot_one function, one palette, one DMR overdraw method.
# ============================================================================
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")


OUT_SMN2      <- "results/smn2_locus_final"
OUT_SENSITIVE <- "results/smn2_sensitive_locus"
OUT_CANDIDATES <- "results/thesis_figures/locus_candidates"
OUT_SMN2_ALL   <- "results/thesis_figures/locus_smn2"
dir.create(OUT_SMN2,      recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_SENSITIVE, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_CANDIDATES, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_SMN2_ALL, recursive=TRUE, showWarnings=FALSE)


COND_COLS <- c(
  ASO_CTRL      = "#1F3A5F",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#C0392B",
  Scramble_VPA  = "#F0A500"
)
DMR_COL  <- "#E31A1C"
WIN_SIZE <- 300
MIN_COV  <- 3


CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond_a="Scramble_CTRL", cond_b="ASO_CTRL",
       label="ASO effect (CTRL background)"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       cond_a="Scramble_CTRL", cond_b="Scramble_VPA",
       label="VPA effect (HDAC inhibitor)"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       cond_a="Scramble_VPA", cond_b="ASO_VPA",
       label="ASO effect on VPA background"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       cond_a="ASO_CTRL", cond_b="ASO_VPA",
       label="VPA effect on ASO background")
)


LOCI_COORDS <- list(
  SMN1 = list(chr="chr5", start=70924941, end=70953015, strand="+"),
  SMN2 = list(chr="chr5", start=70049638, end=70078522, strand="+")
)
EXONS <- list(
  SMN1 = data.frame(
    exon=c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start=c(70925030,70938807,70941357,70942326,70942686,
            70944627,70946033,70951913,70952411),
    end=c(70925158,70938878,70941476,70942526,70942838,
          70944722,70946143,70951966,70952984),
    is_target=c(F,F,F,F,F,F,F,T,F)
  ),
  SMN2 = data.frame(
    exon=c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
    start=c(70049638,70063415,70065965,70066934,70067294,
            70069235,70070641,70076521,70077019),
    end=c(70049766,70063486,70066084,70067134,70067446,
          70069330,70070751,70076574,70077592),
    is_target=c(F,F,F,F,F,F,F,T,F)
  )
)

build_smn_gff <- function() {
  rows <- list()
  for (ln in names(LOCI_COORDS)) {
    l <- LOCI_COORDS[[ln]]
    rows[[length(rows)+1]] <- data.frame(
      chr=l$chr, start=l$start, end=l$end,
      strand=l$strand, type="gene", name=ln)
    ex <- EXONS[[ln]]
    for (i in seq_len(nrow(ex)))
      rows[[length(rows)+1]] <- data.frame(
        chr=l$chr, start=ex$start[i], end=ex$end[i],
        strand=l$strand, type="exon",
        name=sprintf("%s_%s", ln, ex$exon[i]))
  }
  df <- do.call(rbind, rows)
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          strand=df$strand, type=df$type, name=df$name)
}

build_gene_gff <- function(locus) {
  GRanges(seqnames=locus$chr,
          ranges=IRanges(locus$gene_start, locus$gene_end),
          strand=locus$strand, type="gene", name=locus$name)
}


plot_one <- function(meth_a, meth_b, ct, region, gff, title,
                     dmrs=NULL, annotation=NULL,
                     exon_labels=NULL, flank=2000) {
  meth_a <- meth_a[meth_a$readsN >= MIN_COV]
  meth_b <- meth_b[meth_b$readsN >= MIN_COV]

  plotLocalMethylationProfile(
    methylationData1 = meth_a,
    methylationData2 = meth_b,
    region           = region,
    DMRs             = NULL,
    conditionsNames  = c(ct$cond_a, ct$cond_b),
    gff              = gff,
    windowSize       = WIN_SIZE,
    context          = "CG",
    col              = c(COND_COLS[ct$cond_a], COND_COLS[ct$cond_b]),
    main             = title,
    plotMeanLines    = TRUE,
    plotPoints       = TRUE
  )

  # DMR overdraw — red boxes at top of methylation panel
  # usr y-range is -0.256 to 1.256; methylation occupies 0 to 1
  # so ybottom=0.88, ytop=1.00 places boxes at the top of the data area
  if (!is.null(dmrs) && length(dmrs) > 0 && length(dmrs) <= 30) {
    dmr_df <- as.data.frame(dmrs)
    pad <- 500
    for (i in seq_len(nrow(dmr_df))) {
      rect(xleft   = dmr_df$start[i] - pad,
           xright  = dmr_df$end[i]   + pad,
           ybottom = 1.02,
           ytop    = 1.06,
           col     = DMR_COL,
           border  = DMR_COL,
           lwd     = 0.5, xpd=NA)
    }
  }

  # Exon labels for SMN plots
  if (!is.null(exon_labels)) {
    for (i in seq_len(nrow(exon_labels))) {
      mtext(exon_labels$exon[i], side=1,
            at=(exon_labels$start[i]+exon_labels$end[i])/2,
            line=-1.5, cex=0.45,
            col=if(exon_labels$is_target[i]) DMR_COL else "black",
            font=if(exon_labels$is_target[i]) 2 else 1)
    }
  }

  if (!is.null(annotation))
    mtext(annotation, side=3, line=0.2, cex=0.7, col="grey40", adj=1)
}


message("Loading methylation cache...")
meth_pooled <- readRDS("results/dmr/meth_pooled_cache.rds")

message("Loading DMR results...")
dmr_results <- list()
for (ct in CONTRASTS) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (file.exists(rds)) {
    dmr_results[[ct$name]] <- readRDS(rds)
    message("  ", ct$name, ": ", length(dmr_results[[ct$name]]))
  }
}

message("Loading sensitive DMRs...")
SENSITIVE_DMRS <- list()
sens_csv <- "results/smn2_local_dmr/SMN2_sensitive_DMRs_all_contrasts.csv"
if (file.exists(sens_csv)) {
  sdf <- read.csv(sens_csv)
  sdf$contrast <- gsub('"', '', sdf$contrast)
  for (ct in unique(sdf$contrast)) {
    sub <- sdf[sdf$contrast==ct,]
    if (nrow(sub)>0)
      SENSITIVE_DMRS[[ct]] <- GRanges(sub$seqnames,
        IRanges(sub$start, sub$end),
        regionType=sub$regionType, pValue=sub$pValue)
  }
  message("  Loaded ", length(SENSITIVE_DMRS), " contrast sets")
}

SMN_GFF  <- build_smn_gff()
SMN2_REGION <- GRanges("chr5", IRanges(70040000, 70092000))

get_dmrs <- function(ct_name, region) {
  if (is.null(dmr_results[[ct_name]])) return(NULL)
  subsetByOverlaps(dmr_results[[ct_name]], region)
}


# SECTION 1: SMN2 locus (masked, all 4 contrasts)

message("\n=== SECTION 1: SMN2 locus plots ===")
for (ct in CONTRASTS) {
  fname <- paste0("SMN_locus_masked_", ct$name, ".pdf")
  message("  ", fname)
  cairo_pdf(file.path(OUT_SMN2, fname), width=11, height=7,
            bg="white", onefile=TRUE)
  par(mfrow=c(1,1), bg="white", col.axis="black",
      col.lab="black", col.main="black", fg="black")
  dmrs <- get_dmrs(ct$name, SMN2_REGION)
  tryCatch(
    plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
             ct, SMN2_REGION, SMN_GFF,
             title=sprintf("SMN2: %s vs %s\n%s", ct$cond_a, ct$cond_b, ct$label),
             dmrs=dmrs, exon_labels=EXONS[["SMN2"]]),
    error=function(e) { plot.new(); text(0.5,0.5,e$message,col="red") }
  )
  dev.off()
}


# SECTION 2: SMN2 sensitive scan locus (all 4 contrasts)

message("\n=== SECTION 2: SMN2 sensitive scan plots ===")
for (ct in CONTRASTS) {
  fname <- paste0("SMN2_sensitive_", ct$name, ".pdf")
  message("  ", fname)
  cairo_pdf(file.path(OUT_SENSITIVE, fname), width=11, height=7,
            bg="white", onefile=TRUE)
  par(mfrow=c(1,1), bg="white", col.axis="black",
      col.lab="black", col.main="black", fg="black")
  sens_dmrs <- SENSITIVE_DMRS[[ct$name]]
  if (!is.null(sens_dmrs)) sens_dmrs <- subsetByOverlaps(sens_dmrs, SMN2_REGION)
  tryCatch(
    plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
             ct, SMN2_REGION, SMN_GFF,
             title=sprintf("SMN2 sensitive scan: %s vs %s\n%s",
                           ct$cond_a, ct$cond_b, ct$label),
             dmrs=sens_dmrs, exon_labels=EXONS[["SMN2"]],
             annotation="Sensitive DMR parameters: minDiff=5%, p<0.05, minCyto=3"),
    error=function(e) { plot.new(); text(0.5,0.5,e$message,col="red") }
  )
  dev.off()
}


# SECTION 3: Gene locus plots helper function

plot_gene_locus <- function(locus, out_dir) {
  out_path <- file.path(out_dir, paste0(locus$name, "_4contrasts.pdf"))
  message("  ", basename(out_path))
  locus_gr <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  gff <- build_gene_gff(locus)
  cairo_pdf(out_path, width=11, height=6, bg="white", onefile=TRUE)
  for (ct in CONTRASTS) {
    par(bg="white", col.axis="black", col.lab="black",
        col.main="black", fg="black", mar=c(5,4,4,2)+0.1)
    dmrs <- get_dmrs(ct$name, locus_gr)
    tryCatch(
      plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
               ct, locus_gr, gff,
               title=sprintf("%s: %s vs %s (%s)\n%s",
                             locus$name, ct$cond_a, ct$cond_b,
                             ct$label, locus$label),
               dmrs=dmrs, annotation=locus$annotation),
      error=function(e) {
        plot.new()
        text(0.5,0.5,paste0("failed: ",conditionMessage(e)),col="red",cex=0.9)
      }
    )
  }
  dev.off()
}


# SECTION 4: SEMA3C and SIGMAR1

message("\n=== SECTION 4: SEMA3C and SIGMAR1 ===")
example_loci <- list(
  list(name="SEMA3C", chr="chr7", start=80799100, end=80829100,
       gene_start=80765501, gene_end=80815949, strand="+",
       label="SEMA3C — VPA-driven example (axon guidance, chr7)",
       annotation="Intron | hypo | VPA-driven"),
  list(name="SIGMAR1", chr="chr9", start=34625116, end=34655116,
       gene_start=34634739, gene_end=34637844, strand="+",
       label="SIGMAR1 — ASO-driven example (sigma receptor, chr9)",
       annotation="Promoter | hyper | ASO-driven")
)
for (locus in example_loci) plot_gene_locus(locus, OUT_CANDIDATES)


# SECTION 5: Original 7 survivors (pairwise reclassified)

message("\n=== SECTION 5: 7 survivor loci ===")
survivor_loci <- list(
  # ASO context-dependent
  list(name="RELL2", chr="chr5", start=141600122, end=141630122,
       gene_start=141589123, gene_end=141689422, strand="+",
       label="RELL2 — ASO context-dependent (promoter, chr5)",
       annotation="Promoter (1-2kb) | deviation=-0.381 | ASO context-dependent"),
  list(name="DDIT4L", chr="chr4", start=100157963, end=100187963,
       gene_start=100141264, gene_end=100241563, strand="+",
       label="DDIT4L — ASO context-dependent (promoter, chr4)",
       annotation="Promoter (<=1kb) | deviation=-0.275 | ASO context-dependent"),
  list(name="MRPS2", chr="chr9", start=135477616, end=135507616,
       gene_start=135453567, gene_end=135553866, strand="+",
       label="MRPS2 — ASO context-dependent (promoter, chr9)",
       annotation="Promoter (<=1kb) | deviation=-0.252 | ASO context-dependent"),
  list(name="GNG14", chr="chr19", start=12676918, end=12706918,
       gene_start=12641619, gene_end=12742218, strand="+",
       label="GNG14 — ASO context-dependent (intron, chr19)",
       annotation="Intron | deviation=-0.208 | ASO context-dependent"),
  # VPA context-dependent
  list(name="RNA5S13", chr="chr1", start=228662668, end=228692668,
       gene_start=228587769, gene_end=228688068, strand="+",
       label="RNA5S13 — VPA context-dependent (promoter, chr1)",
       annotation="Promoter (<=1kb) | deviation=+0.284 | VPA context-dependent"),
  list(name="KIAA1656", chr="chr22", start=30349684, end=30379684,
       gene_start=30331935, gene_end=30432234, strand="+",
       label="KIAA1656 — VPA context-dependent (intron, chr22)",
       annotation="Intron | deviation=-0.215 | VPA context-dependent"),
  list(name="TCEAL4", chr="chrX", start=103579302, end=103609302,
       gene_start=103536053, gene_end=103636352, strand="+",
       label="TCEAL4 — VPA context-dependent (promoter, chrX)",
       annotation="Promoter (<=1kb) | deviation=+0.244 | VPA context-dependent")
)
for (locus in survivor_loci) plot_gene_locus(locus, OUT_CANDIDATES)


# SECTION 6: New pairwise candidates

message("\n=== SECTION 6: New pairwise candidates ===")
new_candidates <- list(
  # ASO context-dependent
  list(name="IRF8", chr="chr16", start=85879034, end=85909034,
       gene_start=85865935, gene_end=85966234, strand="+",
       label="IRF8 — ASO context-dependent (transcription factor, chr16)",
       annotation="Promoter | score=9 | neuroinflammation regulator"),
  list(name="USP27X", chr="chrX", start=49909302, end=49939302,
       gene_start=49832753, gene_end=49933052, strand="+",
       label="USP27X — ASO context-dependent (deubiquitinase, chrX)",
       annotation="Downstream | score=7 | ASO context-dependent"),
  list(name="USP7", chr="chr16", start=8912084, end=8942084,
       gene_start=8840635, gene_end=8940934, strand="+",
       label="USP7 — ASO context-dependent (deubiquitinase/chromatin, chr16)",
       annotation="Distal intergenic | score=6 | ASO context-dependent"),
  list(name="KDM1A", chr="chr1", start=23034568, end=23064568,
       gene_start=23035569, gene_end=23135868, strand="+",
       label="KDM1A — ASO context-dependent (histone demethylase LSD1, chr1)",
       annotation="3 UTR | score=4 | histone demethylase"),
  list(name="PAX5_ASO", chr="chr9", start=37046116, end=37076116,
       gene_start=36992067, gene_end=37092366, strand="+",
       label="PAX5 — ASO context-dependent (transcription factor, chr9)",
       annotation="Intron | score=6 | ASO context-dependent"),
  # VPA context-dependent
  list(name="PAX5_VPA", chr="chr9", start=37046116, end=37076116,
       gene_start=36992067, gene_end=37092366, strand="+",
       label="PAX5 — VPA context-dependent (transcription factor, chr9)",
       annotation="Intron | score=7 | VPA context-dependent"),
  list(name="CAMK2A", chr="chr5", start=150245072, end=150275072,
       gene_start=150164623, gene_end=150264922, strand="+",
       label="CAMK2A — VPA context-dependent (calcium kinase, chr5)",
       annotation="Intron | score=7 | synaptic/motor neuron"),
  list(name="EPHB1", chr="chr3", start=135030969, end=135060969,
       gene_start=135041420, gene_end=135141719, strand="+",
       label="EPHB1 — VPA context-dependent (axon guidance, chr3)",
       annotation="Intron | score=7 | axon guidance"),
  list(name="ZDHHC22", chr="chr14", start=77081022, end=77111022,
       gene_start=77092523, gene_end=77192822, strand="+",
       label="ZDHHC22 — VPA context-dependent (palmitoyl transferase, chr14)",
       annotation="Promoter | score=6 | VPA context-dependent")
)
for (locus in new_candidates) plot_gene_locus(locus, OUT_CANDIDATES)

message("\n=== ALL DONE ===")
message("SMN2 locus:      ", OUT_SMN2)
message("SMN2 sensitive:  ", OUT_SENSITIVE)
message("All gene loci:   ", OUT_CANDIDATES)
