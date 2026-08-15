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

# Coordinates centred on actual DMR ± 5kb
LOCI <- list(
  list(name="HUWE1",  chr="chrX",  start=53662000, end=53763000,
       desc="E3 ubiquitin ligase — neural development | ASO=29% ASO+VPA=53%"),
  list(name="GFRA2",  chr="chr8",  start=21763000, end=21864000,
       desc="GDNF receptor alpha 2 — motor neuron survival | ASO=42% ASO+VPA=58%"),
  list(name="ROCK1",  chr="chr18", start=19675000, end=19776000,
       desc="Rho kinase 1 — motor neuron survival | ASO=29% ASO+VPA=46%"),
  list(name="TSHZ1",  chr="chr18", start=75254000, end=75355000,
       desc="Teashirt zinc finger homeobox 1 — neural transcription factor | ASO=21% ASO+VPA=46%"),
  list(name="CD38",   chr="chr4",  start=15792000, end=15894000,
       desc="CD38 — NAD+ metabolism, neurodegeneration | ASO=21% ASO+VPA=46%"),
  list(name="SLC32A1",chr="chr20", start=38670000, end=38771000,
       desc="GABA vesicular transporter — inhibitory synapse | ASO=46% ASO+VPA=57%")
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

for (locus in LOCI) {
  message("Plotting: ", locus$name, " (", locus$chr, ")")
  pooled <- lapply(CONDITIONS, load_cpg, chr=locus$chr)
  names(pooled) <- CONDITIONS
  pooled <- Filter(Negate(is.null), pooled)

  region <- GRanges(locus$chr, IRanges(locus$start, locus$end))

  df_all <- do.call(rbind, lapply(names(pooled), function(cond) {
    prof <- computeMethylationProfile(pooled[[cond]], region, 300, "CG")
    df <- as.data.frame(prof)
    df$meth <- df$sumReadsM / df$sumReadsN
    df$pos  <- (df$start + df$end) / 2
    df$cond <- cond
    df[!is.na(df$meth) & df$sumReadsN >= 3, ]
  }))
  df_all$cond <- factor(df_all$cond, levels=CONDITIONS)

  df_all$linetype <- ifelse(grepl("Scramble", df_all$cond), "dashed", "solid")
  p <- ggplot(df_all, aes(x=pos/1e6, y=meth, colour=cond, linetype=cond)) +
    geom_line(linewidth=1.1) +
    scale_colour_manual(values=COND_COLOURS, name=NULL) +
    scale_linetype_manual(values=c(ASO_CTRL="solid", Scramble_CTRL="dashed",
                                   ASO_VPA="solid", Scramble_VPA="dashed"), name=NULL) +
    scale_y_continuous(limits=c(0,1), labels=scales::percent_format(1)) +
    theme_classic(base_size=13) +
    theme(legend.position="right",
          plot.title=element_text(face="bold", size=12),
          plot.subtitle=element_text(colour="grey30", size=9)) +
    labs(title=paste0(locus$name, " — ", locus$desc),
         subtitle=paste0(locus$chr,":",format(locus$start,big.mark=","),
                         "-",format(locus$end,big.mark=","),
                         " | 300bp bins | gap = methylation difference between conditions"),
         x=paste0(locus$chr," position (Mb)"),
         y="CpG methylation (proportion)")

  fname <- file.path(OUT_DIR, paste0("lowres_dip_",locus$name,"_4cond.pdf"))
  ggsave(fname, p, width=11, height=5, device=cairo_pdf)
  message("saved: ", basename(fname))
}
message("Done.")
