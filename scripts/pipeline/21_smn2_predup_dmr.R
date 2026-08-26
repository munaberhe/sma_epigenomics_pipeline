#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT     <- "results/smn2_predup_comparison"
CX_DIR  <- "results/smn2_predup_cx"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

SMN2_REG <- GRanges("chr5", IRanges(70040000, 70090000))

CONDITIONS <- c("ASO_CTRL","ASO_VPA","Scramble_CTRL","Scramble_VPA")
COND_COLS  <- c(ASO_CTRL="#1F3A5F", Scramble_CTRL="#6B7280",
                ASO_VPA="#C0392B",  Scramble_VPA="#F0A500")

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",
       cond1="ASO_CTRL", cond2="Scramble_CTRL"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",
       cond1="Scramble_VPA", cond2="Scramble_CTRL"),
  list(name="ASO_VPA_vs_Scramble_VPA",
       cond1="ASO_VPA", cond2="Scramble_VPA"),
  list(name="ASO_VPA_vs_ASO_CTRL",
       cond1="ASO_VPA", cond2="ASO_CTRL")
)

message("Loading pre-extracted SMN2 CX data...")
load_smn2 <- function(condition) {
  message("  ", condition)
  grs <- lapply(1:3, function(r) {
    # match filename pattern
    pattern <- paste0(condition, "_", r, "_*_smn2.txt")
    f <- list.files(CX_DIR, pattern=glob2rx(pattern), full.names=TRUE)
    if (length(f) == 0) {
      message("    missing: ", pattern)
      return(NULL)
    }
    d <- read.table(f[1], header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","countM","countU","context","tri"),
      colClasses=c("character","integer","character","integer",
                   "integer","character","character"))
    if (nrow(d) == 0) return(NULL)
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos, d$pos),
            strand=d$strand, readsM=d$countM,
            readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  grs <- Filter(Negate(is.null), grs)
  if (length(grs) == 0) return(NULL)
  poolMethylationDatasets(GRangesList(grs))
}

meth_predup <- lapply(CONDITIONS, load_smn2)
names(meth_predup) <- CONDITIONS

# Coverage report
for (cond in CONDITIONS) {
  m <- meth_predup[[cond]]
  if (is.null(m)) { message(cond, ": NO DATA"); next }
  message(cond, ": ", length(m), " CpGs | median cov: ",
          round(median(m$readsN),1), " | >=5x: ", sum(m$readsN >= 5),
          " | >=10x: ", sum(m$readsN >= 10))
}

message("\nStage 1: canonical DMR parameters...")
dmr_predup <- list()
for (ct in CONTRASTS) {
  m1 <- meth_predup[[ct$cond1]]
  m2 <- meth_predup[[ct$cond2]]
  if (is.null(m1) || is.null(m2)) next
  tryCatch({
    dmrs <- computeDMRs(m1, m2,
      method="bins", context="CG",
      binSize=300, minCytosinesCount=4,
      minProportionDifference=0.20,
      minReadsPerCytosine=4, pValueThreshold=0.01,
      minGap=300, test="score")
    dmr_predup[[ct$name]] <- dmrs
    message("  ", ct$name, ": ", length(dmrs), " DMRs")
  }, error=function(e) message("  ", ct$name, ": ERROR - ", e$message))
}
saveRDS(dmr_predup, file.path(OUT, "dmr_canonical_predup_SMN2.rds"))

message("\nStage 2: sensitive sweep...")
PARAM_SETS <- list(
  list(label="5pct_100bp", minDiff=0.05, binSize=100, minCyto=3, pval=0.05),
  list(label="2pct_50bp",  minDiff=0.02, binSize=50,  minCyto=2, pval=0.05),
  list(label="1pct_50bp",  minDiff=0.01, binSize=50,  minCyto=2, pval=0.10),
  list(label="5pct_200bp", minDiff=0.05, binSize=200, minCyto=3, pval=0.05),
  list(label="5pct_500bp", minDiff=0.05, binSize=500, minCyto=3, pval=0.05)
)
dmr_sensitive <- list()
for (ct in CONTRASTS) {
  m1 <- meth_predup[[ct$cond1]]
  m2 <- meth_predup[[ct$cond2]]
  if (is.null(m1) || is.null(m2)) next
  for (ps in PARAM_SETS) {
    key <- paste0(ct$name, "_", ps$label)
    tryCatch({
      dmrs <- computeDMRs(m1, m2,
        regions=SMN2_REG, context="CG", method="bins",
        binSize=ps$binSize, minCytosinesCount=ps$minCyto,
        minProportionDifference=ps$minDiff,
        minReadsPerCytosine=3, pValueThreshold=ps$pval,
        test="score")
      dmr_sensitive[[key]] <- dmrs
      message("  ", key, ": ", length(dmrs), " DMRs")
    }, error=function(e) message("  ", key, ": ERROR - ", e$message))
  }
}
saveRDS(dmr_sensitive, file.path(OUT, "dmr_sensitive_predup_SMN2.rds"))

message("\nPlotting comparison...")
dedup_cache <- readRDS("results/dmr/meth_pooled_cache.rds")
dedup_smn2  <- lapply(dedup_cache, function(m) subsetByOverlaps(m, SMN2_REG))

make_profile <- function(meth_list, title, min_cov=3) {
  df <- do.call(rbind, lapply(names(meth_list), function(cond) {
    m <- meth_list[[cond]]
    if (is.null(m)) return(NULL)
    m <- m[!is.na(m$readsN) & m$readsN >= min_cov]
    if (length(m) == 0) return(NULL)
    data.frame(pos=start(m), meth=m$readsM/m$readsN,
               coverage=m$readsN, condition=cond)
  }))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  ggplot(df, aes(x=pos, y=meth, colour=condition)) +
    geom_point(size=0.8, alpha=0.5) +
    geom_smooth(method="loess", span=0.5, se=FALSE, linewidth=1.2) +
    scale_colour_manual(values=COND_COLS) +
    scale_y_continuous(limits=c(0,1), labels=scales::percent_format(1)) +
    scale_x_continuous(labels=function(x) sprintf("%.3fMb", x/1e6)) +
    geom_vline(xintercept=c(70076521,70076574),
               linetype="dashed", colour="#E31A1C", linewidth=0.5) +
    annotate("text", x=70076547, y=0.05, label="E7",
             colour="#E31A1C", size=3, fontface="bold") +
    labs(title=title, x="chr5 position", y="CpG methylation", colour=NULL) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"), legend.position="right",
          panel.grid.major.y=element_line(colour="grey92"))
}

p_dedup  <- make_profile(dedup_smn2,  "WITH deduplication (canonical)", min_cov=3)
p_predup <- make_profile(meth_predup, "WITHOUT deduplication (pre-dedup)", min_cov=3)

if (!is.null(p_dedup) && !is.null(p_predup)) {
  combined <- (p_dedup / p_predup) +
    plot_annotation(
      title="SMN2 coverage gap test: dedup vs pre-dedup",
      caption=paste0("Red dashed lines = Exon 7 boundaries.\n",
        "If gap closes in pre-dedup panel: deduplication caused the coverage gap.\n",
        "If gap persists: sequencing capture limitation."),
      theme=theme(plot.title=element_text(face="bold", size=13),
                  plot.caption=element_text(colour="grey50", size=9))
    )
  ggsave(file.path(OUT, "SMN2_predup_vs_dedup_comparison.pdf"),
         combined, width=12, height=10, device=cairo_pdf)
  ggsave(file.path(OUT, "SMN2_predup_vs_dedup_comparison.png"),
         combined, width=12, height=10, dpi=150)
  message("Saved comparison plot")
}
message("\nDone. Outputs in: ", OUT)
