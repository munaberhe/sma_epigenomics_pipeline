.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(DMRcaller)
  library(ggplot2)
  library(patchwork)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/lowres_profiles"

COND_COLOURS <- c(ASO_CTRL="#1F3A5F", Scramble_CTRL="#6B7280",
                  ASO_VPA="#C0392B", Scramble_VPA="#D4A017")
CONDITIONS <- c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA")
BY_CHR <- "results/alignments/bs/by_chr"

LOCI <- list(
  list(name="chrX_overview", chr="chrX",
       start=1, end=156040895, winsize=500000,
       desc="Full chrX — 500kb bins | diffuse ASO hypomethylation"),
  list(name="chrX_hotspot", chr="chrX",
       start=5000000, end=55000000, winsize=50000,
       desc="chrX 5-55Mb hotspot — 50kb bins | peak DMR density regions"),
  list(name="SYP",  chr="chrX", start=49091000, end=49291000, winsize=300,
       desc="Synaptophysin — presynaptic vesicle protein"),
  list(name="EFNB1",chr="chrX", start=68790000, end=68990000, winsize=300,
       desc="Ephrin B1 — axon guidance / neural circuits")
)

load_cpg <- function(condition, chr) {
  files <- file.path(BY_CHR,
    sprintf("%s_%d_%s.CpG_report.txt.gz", condition, 1:3, chr))
  files <- files[file.exists(files)]
  if (length(files)==0) return(NULL)
  grs <- lapply(files, function(f) {
    d <- read.table(gzfile(f), header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","M","U","ctx","tri"))
    d <- d[d$ctx=="CG",]
    GRanges(d$chr, IRanges(d$pos,d$pos), strand=d$strand,
            readsM=d$M, readsN=d$M+d$U, context=d$ctx,
            trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

message("Loading chrX CpG data...")
pooled <- lapply(CONDITIONS, load_cpg, chr="chrX")
names(pooled) <- CONDITIONS
pooled <- Filter(Negate(is.null), pooled)

for (locus in LOCI) {
  message("Plotting: ", locus$name)
  region <- GRanges(locus$chr, IRanges(locus$start, locus$end))

  df_all <- do.call(rbind, lapply(names(pooled), function(cond) {
    prof <- computeMethylationProfile(pooled[[cond]], region, locus$winsize, "CG")
    df <- as.data.frame(prof)
    df$meth <- df$sumReadsM / df$sumReadsN
    df$pos  <- (df$start + df$end) / 2
    df$cond <- cond
    df[!is.na(df$meth) & df$sumReadsN >= 3, ]
  }))
  df_all$cond <- factor(df_all$cond, levels=CONDITIONS)

  p <- ggplot(df_all, aes(x=pos/1e6, y=meth, colour=cond, linetype=cond)) +
    geom_line(linewidth=1.0) +
    scale_colour_manual(values=COND_COLOURS, name=NULL) +
    scale_linetype_manual(
      values=c(ASO_CTRL="solid", Scramble_CTRL="dashed",
               ASO_VPA="solid", Scramble_VPA="dashed"), name=NULL) +
    scale_y_continuous(limits=c(0,1),
                       
                       breaks=seq(0,1,0.2)) +
    scale_x_continuous(labels=function(x) paste0(round(x,1)," Mb")) +
    theme_classic(base_size=13) +
    theme(legend.position="right",
          plot.title=element_text(face="bold", size=12),
          plot.subtitle=element_text(colour="grey30", size=9),
          panel.grid.major.y=element_line(colour="grey92")) +
    labs(title=paste0(locus$name," — ",locus$desc),
         subtitle=paste0("chrX | solid=ASO conditions, dashed=Scramble controls | gap = methylation difference"),
         x="chrX position (Mb)", y="CpG methylation (proportion)")

  fname <- file.path(OUT_DIR, paste0("lowres_chrX_",locus$name,"_300bp.pdf"))
  ggsave(fname, p, width=12, height=5, device=cairo_pdf)
  message("saved: ", basename(fname))
}
message("Done.")
