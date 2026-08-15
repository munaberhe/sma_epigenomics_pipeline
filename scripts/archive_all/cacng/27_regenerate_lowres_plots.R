.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges); library(DMRcaller)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
source("scripts/cacng/26_lowres_helper.R")

OUT_DIR <- "results/lowres_profiles"
BY_CHR  <- "results/alignments/bs/by_chr"
CONDITIONS <- c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA")

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

# ---- Define all loci ----
LOCI <- list(
  list(name="HUWE1", chr="chrX", start=53638253, end=53678852,
       genes=data.frame(name=c("HUWE1"), start=c(53658253), end=c(53658852), focal=c(TRUE)),
       title="HUWE1: E3 ubiquitin ligase — neural development",
       subtitle="chrX | ASO=29% hypo, ASO+VPA=53% hypo | solid=ASO, dashed=Scramble"),
  list(name="GFRA2", chr="chr8", start=21733107, end=21773706,
       genes=data.frame(name=c("GFRA2"), start=c(21753107), end=c(21753706), focal=c(TRUE)),
       title="GFRA2: GDNF receptor alpha 2 — motor neuron survival",
       subtitle="chr8 | ASO=42% hypo, ASO+VPA=58% hypo | solid=ASO, dashed=Scramble"),
  list(name="ROCK1", chr="chr18", start=20953790, end=20994089,
       genes=data.frame(name=c("ROCK1"), start=c(20973790), end=c(20974089), focal=c(TRUE)),
       title="ROCK1: Rho kinase 1 — motor neuron survival",
       subtitle="chr18 | ASO=29% hypo, ASO+VPA=46% hypo | solid=ASO, dashed=Scramble"),
  list(name="TSHZ1", chr="chr18", start=75194390, end=75235289,
       genes=data.frame(name=c("TSHZ1"), start=c(75214390), end=c(75215289), focal=c(TRUE)),
       title="TSHZ1: Teashirt zinc finger — neural transcription factor",
       subtitle="chr18 | ASO=21% hypo, ASO+VPA=46% hypo | solid=ASO, dashed=Scramble"),
  list(name="CD38", chr="chr4", start=15767464, end=15807763,
       genes=data.frame(name=c("CD38"), start=c(15787464), end=c(15787763), focal=c(TRUE)),
       title="CD38: NAD+ metabolism — neurodegeneration link",
       subtitle="chr4 | ASO=21% hypo, ASO+VPA=46% hypo | solid=ASO, dashed=Scramble"),
  list(name="SLC32A1", chr="chr20", start=38706409, end=38746708,
       genes=data.frame(name=c("SLC32A1"), start=c(38726409), end=c(38726708), focal=c(TRUE)),
       title="SLC32A1: GABA vesicular transporter — inhibitory synapse",
       subtitle="chr20 | ASO=46% hypo, ASO+VPA=57% hypo | solid=ASO, dashed=Scramble"),
  list(name="CACNG_cluster", chr="chr17", start=66909154, end=67081753, dmr_start=66911154, dmr_end=67079753,
       genes=data.frame(
         name=c("CACNG5","CACNG4","CACNG1"),
         start=c(66835117, 66964707, 67044554),
         end=c(66894751, 67033398, 67056797),
         focal=c(TRUE, TRUE, TRUE)),
       title="chr17 CACNG cluster: VPA reverses ASO hypermethylation",
       subtitle="CACNG5, CACNG4, CACNG1 — AMPA receptor subunits | chr17:66.8-67.1 Mb | solid=ASO, dashed=Scramble"),

  list(name="CACNG1", chr="chr17", start=67074454, end=67084753,
       dmr_start=67079454, dmr_end=67079753,
       genes=data.frame(name=c("CACNG1 DMR"), start=c(67079454), end=c(67079753), focal=c(TRUE)),
       title="CACNG1: AMPA receptor subunit — combination-exclusive",
       subtitle="chr17 | no ASO effect alone | ASO+VPA=+28% hypo"),
  list(name="CACNG4", chr="chr17", start=66992554, end=67003153,
       dmr_start=66997554, dmr_end=66998153,
       genes=data.frame(name=c("CACNG4 DMR"), start=c(66997554), end=c(66998153), focal=c(TRUE)),
       title="CACNG4: AMPA receptor subunit — combination-exclusive",
       subtitle="chr17 | no ASO effect alone | ASO+VPA=+26% hypo"),
  list(name="CACNG5", chr="chr17", start=66906154, end=66916753,
       dmr_start=66911154, dmr_end=66911753,
       genes=data.frame(name=c("CACNG5 DMR"), start=c(66911154), end=c(66911753), focal=c(TRUE)),
       title="CACNG5: AMPA receptor subunit — VPA reverses ASO hypermethylation",
       subtitle="chr17 | ASO=-24% hyper | ASO+VPA=+25% hypo"),
  list(name="PRKCA", chr="chr17", start=66306454, end=66316753,
       dmr_start=66311454, dmr_end=66311753,
       genes=data.frame(name=c("PRKCA DMR"), start=c(66311454), end=c(66311753), focal=c(TRUE)),
       title="PRKCA: Protein kinase C alpha — combination-exclusive",
       subtitle="chr17 | no ASO effect alone | ASO+VPA=+27% hypo"),
  list(name="SEMA3C", chr="chr7", start=80808801, end=80819400,
       dmr_start=80813801, dmr_end=80814400,
       genes=data.frame(name=c("SEMA3C DMR"), start=c(80813801), end=c(80814400), focal=c(TRUE)),
       title="SEMA3C: Semaphorin 3C — axon guidance",
       subtitle="chr7 | ASO+VPA strongest hypo hit p=3.6e-69"),
  list(name="MTA1", chr="chr14", start=105423023, end=105433922,
       dmr_start=105428023, dmr_end=105428922,
       genes=data.frame(name=c("MTA1 DMR"), start=c(105428023), end=c(105428922), focal=c(TRUE)),
       title="MTA1: Metastasis associated 1 — chromatin remodelling",
       subtitle="chr14 | ASO+VPA hypo p=2.1e-09"),
  list(name="SYP", chr="chrX", start=49184453, end=49224752,
       genes=data.frame(name=c("SYP"), start=c(49204453), end=c(49204752), focal=c(TRUE)),
       title="SYP: Synaptophysin — presynaptic vesicle protein",
       subtitle="chrX | ASO-specific hypomethylation | solid=ASO, dashed=Scramble"),
  list(name="EFNB1", chr="chrX", start=68838053, end=68878352,
       genes=data.frame(name=c("EFNB1"), start=c(68858053), end=c(68858352), focal=c(TRUE)),
       title="EFNB1: Ephrin B1 — axon guidance",
       subtitle="chrX | ASO-specific hypomethylation | solid=ASO, dashed=Scramble")
)

# load and plot each locus
for (locus in LOCI) {
  message("\n=== ", locus$name, " ===")
  pooled <- lapply(CONDITIONS, load_cpg, chr=locus$chr)
  names(pooled) <- CONDITIONS
  pooled <- Filter(Negate(is.null), pooled)
  if (length(pooled)==0) { message("no data"); next }

  region <- GRanges(locus$chr, IRanges(locus$start, locus$end))
  out_path <- file.path(OUT_DIR, paste0("lowres_v2_",locus$name,"_300bp.pdf"))

  plot_locus_baseR(pooled=pooled, CONDITIONS=CONDITIONS, region=region,
             winsize=300, genes_df=locus$genes,
             dmr_start=locus$dmr_start, dmr_end=locus$dmr_end,
             title=locus$title, out_path=out_path)
}
message("\nAll done. Outputs in: ", OUT_DIR)
