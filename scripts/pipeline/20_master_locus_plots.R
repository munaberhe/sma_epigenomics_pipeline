#!/usr/bin/env Rscript
# ============================================================================
# 20_master_locus_plots.R
# Generates all CpG methylation locus plots for the thesis.
# Sections:
#   1. SMN2 locus (masked reference, 4 pairwise contrasts)
#   2. SMN2 sensitive scan (relaxed parameters, 4 contrasts)
#   3. Gene locus helper function
#   4. Pairwise synergy candidates (8 genes in both ASO and VPA context-dependent lists)
#   5. Pairwise context-dependent candidates (9 genes, ASO or VPA specific)
# All candidate coordinates taken from pairwise_context_scan scored CSVs.
# Usage: sbatch --mem=64G --wrap="Rscript scripts/pipeline/20_master_locus_plots.R"
# ============================================================================
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")


OUT_SMN2      <- "results/smn2_locus_final"
OUT_SENSITIVE <- "results/smn2_sensitive_locus"
OUT_SYNERGY   <- "results/thesis_figures/locus_candidates/01_synergy"
OUT_ASO_RESTR <- "results/thesis_figures/locus_candidates/02_ASO_restricted"
OUT_VPA_RESTR <- "results/thesis_figures/locus_candidates/03_VPA_restricted"
OUT_SMN2_ALL  <- "results/thesis_figures/locus_smn2"
dir.create(OUT_SMN2,      recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_SENSITIVE, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_SYNERGY,   recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_ASO_RESTR, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_VPA_RESTR, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_SMN2_ALL,  recursive=TRUE, showWarnings=FALSE)


COND_COLS <- c(
  ASO_CTRL      = "#1F3A5F",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#C0392B",
  Scramble_VPA  = "#F0A500"
)
DMR_COL  <- "#C0392B"
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
    # label the DMR track so readers know what the boxes are
    text(x=par("usr")[1], y=1.04, labels="DMR", col=DMR_COL, font=2, cex=0.6, xpd=NA, adj=c(0,0.5))
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

message("Plotting SMN2 locus (masked)...")
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
             title=sprintf("SMN2: %s vs %s\n%s", ct$cond_b, ct$cond_a, ct$label),
             dmrs=dmrs, exon_labels=EXONS[["SMN2"]]),
    error=function(e) { plot.new(); text(0.5,0.5,e$message,col="red") }
  )
  dev.off()
}


# SECTION 2: SMN2 sensitive scan locus (all 4 contrasts)

message("Plotting SMN2 sensitive scan...")
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
    # compute per-contrast annotation from actual DMR stats
    # locus$annotation contains only the genomic feature (e.g. "Promoter")
    # meth_diff and p are computed fresh per contrast from the DMR data
    if (!is.null(dmrs) && length(dmrs) > 0) {
      dmr_df <- as.data.frame(dmrs)
      dmr_df$meth_diff <- dmr_df$proportion1 - dmr_df$proportion2
      best <- dmr_df[which.min(dmr_df$pValue), ]
      feature <- gsub(" \\| meth_diff.*", "", locus$annotation)
      ct_annotation <- sprintf("%s | meth_diff=%.3f | p=%.2e | n=%d DMRs in window",
                                feature,
                                best$meth_diff, best$pValue, nrow(dmr_df))
    } else {
      feature <- gsub(" \\| meth_diff.*", "", locus$annotation)
      ct_annotation <- sprintf("%s | no DMR in this contrast", feature)
    }
    tryCatch(
      plot_one(meth_pooled[[ct$cond_a]], meth_pooled[[ct$cond_b]],
               ct, locus_gr, gff,
               title=sprintf("%s: %s", locus$name, ct$label),
               dmrs=dmrs, annotation=ct_annotation),
      error=function(e) {
        plot.new()
        text(0.5,0.5,paste0("failed: ",conditionMessage(e)),col="red",cex=0.9)
      }
    )
  }
  dev.off()
}


# synergy candidates: DMR in both ASO_in_VPA and VPA_in_ASO,
# absent from corresponding single-drug contrast.
# coordinates from results/pairwise_context_scan/all_candidate_coords.csv
# synergy: absent from ASO_alone and VPA_alone, present in both combination contrasts
# verified by check_all_candidates2.R
message("Plotting synergy candidates...")
OUT_CANDIDATES <- OUT_SYNERGY
synergy_loci <- list(
  list(name="GMPR2", chr="chr14",
       start=24219672, end=24249672,
       gene_start=24235374, gene_end=24237152, strand="+",
       label="GMPR2 (synergy, promoter, chr14)",
       annotation="Promoter | meth_diff=-0.238 | p=8.79e-04"),
  list(name="C1QA", chr="chr1",
       start=22622218, end=22652218,
       gene_start=22636762, gene_end=22639608, strand="+",
       label="C1QA (synergy, promoter, chr1)",
       annotation="Promoter | meth_diff=-0.262 | p=0.00166"),
  list(name="DDIT4L", chr="chr4",
       start=100176264, end=100206563,
       gene_start=100187875, gene_end=100191563, strand="-",
       label="DDIT4L (synergy, promoter, chr4)",
       annotation="Promoter | meth_diff=-0.201 | p=5.36e-04"),
  list(name="IL13", chr="chr5",
       start=132643172, end=132673172,
       gene_start=132658173, gene_end=132660332, strand="+",
       label="IL13 (synergy, promoter, chr5)",
       annotation="Promoter | meth_diff=-0.250 | p=8.25e-04"),
  list(name="RELL2", chr="chr5",
       start=141624272, end=141654272,
       gene_start=141637442, gene_end=141641053, strand="+",
       label="RELL2 (synergy, promoter, chr5)",
       annotation="Promoter | meth_diff=-0.325 | p=4.04e-04"),
  list(name="IRF8", chr="chr16",
       start=85901084, end=85931084,
       gene_start=85914043, gene_end=85921371, strand="+",
       label="IRF8 (synergy, promoter, chr16)",
       annotation="Promoter | meth_diff=-0.545 | p=7.26e-04"),
  list(name="HDAC4", chr="chr2",
       start=238987750, end=239017750,
       gene_start=239053046, gene_end=239066981, strand="-",
       label="HDAC4 (synergy, exon, chr2)",
       annotation="Exon | meth_diff=-0.321 | p=7.88e-04"),
  list(name="BDNF", chr="chr11",
       start=27705226, end=27735226,
       gene_start=27654893, gene_end=27719747, strand="-",
       label="BDNF (synergy, promoter, chr11)",
       annotation="Promoter | meth_diff=-0.269 | p=1.60e-04"),
  list(name="DLG4", chr="chr17",
       start=7176204, end=7206204,
       gene_start=7190612, gene_end=7192210, strand="-",
       label="DLG4 (synergy, promoter, chr17)",
       annotation="Promoter | meth_diff=-0.253 | p=1.56e-05"),
  list(name="KDM1A", chr="chr1",
       start=23070718, end=23100718,
       gene_start=23079363, gene_end=23083638, strand="+",
       label="KDM1A (synergy, 3UTR, chr1)",
       annotation="3 UTR | meth_diff=-0.468 | p=3.90e-05")
)
for (locus in synergy_loci) plot_gene_locus(locus, OUT_CANDIDATES)

# ASO-background-specific: DMR in ASO_in_VPA only,
# absent from ASO_alone, VPA_alone, and VPA_in_ASO
message("Plotting ASO-background-specific candidates...")
OUT_CANDIDATES <- OUT_ASO_RESTR
aso_bg_loci <- list(
  list(name="TMEM179B", chr="chr11",
       start=62774026, end=62804026,
       gene_start=62789060, gene_end=62790110, strand="+",
       label="TMEM179B (ASO-background-specific, promoter, chr11)",
       annotation="Promoter | meth_diff=-0.439 | p=0.00257"),
  list(name="SRM", chr="chr1",
       start=11037418, end=11067418,
       gene_start=11054584, gene_end=11055204, strand="-",
       label="SRM (ASO-background-specific, intergenic, chr1)",
       annotation="Distal intergenic | meth_diff=0.244 | p=0.00223"),
  list(name="ANKZF1", chr="chr2",
       start=219222850, end=219252850,
       gene_start=219235645, gene_end=219236651, strand="+",
       label="ANKZF1 (ASO-background-specific, exon, chr2)",
       annotation="Exon | meth_diff=-0.365 | p=6.41e-04"),
  list(name="AMT", chr="chr3",
       start=49407070, end=49437070,
       gene_start=49416794, gene_end=49422172, strand="-",
       label="AMT (ASO-background-specific, promoter, chr3)",
       annotation="Promoter | meth_diff=-0.232 | p=2.12e-04"),
  list(name="TAF11L11", chr="chr5",
       start=17587172, end=17617172,
       gene_start=17604177, gene_end=17605377, strand="+",
       label="TAF11L11 (ASO-background-specific, promoter, chr5)",
       annotation="Promoter | meth_diff=0.272 | p=0.00364")
)
for (locus in aso_bg_loci) plot_gene_locus(locus, OUT_CANDIDATES)

# VPA-restricted candidates: DMR in VPA_in_ASO only,
# absent from VPA_alone and ASO_in_VPA
message("Plotting VPA-restricted candidates...")
OUT_CANDIDATES <- OUT_VPA_RESTR
vpa_restricted_loci <- list(
  list(name="AGRN", chr="chr1",
       start=1003018, end=1033018,
       gene_start=1020120, gene_end=1056116, strand="+",
       label="AGRN (VPA-restricted, promoter, chr1)",
       annotation="Promoter | meth_diff=-0.410 | p=1.05e-04"),
  list(name="APH1A", chr="chr1",
       start=150252118, end=150282118,
       gene_start=150267097, gene_end=150267696, strand="-",
       label="APH1A (VPA-restricted, promoter, chr1)",
       annotation="Promoter | meth_diff=-0.440 | p=2.74e-04"),
  list(name="SMN2", chr="chr5",
       start=70073372, end=70103372,
       gene_start=70070664, gene_end=70078522, strand="+",
       label="SMN2 (VPA-restricted, intergenic, chr5)",
       annotation="Distal intergenic | meth_diff=-0.354 | p=2.31e-05"),
  list(name="IRF9", chr="chr14",
       start=24147372, end=24177372,
       gene_start=24162399, gene_end=24166370, strand="+",
       label="IRF9 (VPA-restricted, promoter, chr14)",
       annotation="Promoter | meth_diff=-0.476 | p=3.46e-10"),
  list(name="KANSL1", chr="chr17",
       start=46129404, end=46159404,
       gene_start=46050650, gene_end=46121338, strand="-",
       label="KANSL1 (VPA-restricted, intron, chr17)",
       annotation="Intron | meth_diff=-0.463 | p=0.00186"),
  list(name="KRI1", chr="chr19",
       start=10549768, end=10579768,
       gene_start=10553576, gene_end=10565971, strand="-",
       label="KRI1 (VPA-restricted, promoter, chr19)",
       annotation="Promoter | meth_diff=-0.429 | p=7.49e-04"),
  list(name="MAGEH1", chr="chrX",
       start=55442502, end=55472502,
       gene_start=55452127, gene_end=55453566, strand="+",
       label="MAGEH1 (VPA-restricted, intergenic, chrX)",
       annotation="Distal intergenic | meth_diff=0.430 | p=0.00286")
)
for (locus in vpa_restricted_loci) plot_gene_locus(locus, OUT_CANDIDATES)

message("Done.")
message("SMN2 locus:      ", OUT_SMN2)
message("SMN2 sensitive:  ", OUT_SENSITIVE)
message("Candidates:      ", OUT_CANDIDATES)
