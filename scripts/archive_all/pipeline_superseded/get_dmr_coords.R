.libPaths("~/R/library")
library(GenomicRanges)
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

genes <- list(
  list(name="SEMA3C",   chr="chr7",  start=80760801,  end=80867400),
  list(name="SIGMAR1",  chr="chr9",  start=34589667,  end=34692366),
  list(name="RELL2",    chr="chr5",  start=141589123, end=141689422),
  list(name="DDIT4L",   chr="chr4",  start=100141264, end=100241563),
  list(name="MRPS2",    chr="chr9",  start=135453567, end=135553866),
  list(name="GNG14",    chr="chr19", start=12641619,  end=12742218),
  list(name="RNA5S13",  chr="chr1",  start=228587769, end=228688068),
  list(name="KIAA1656", chr="chr22", start=30331935,  end=30432234),
  list(name="TCEAL4",   chr="chrX",  start=103536053, end=103636352),
  list(name="IRF8",     chr="chr16", start=85865935,  end=85966234),
  list(name="USP27X",   chr="chrX",  start=49832753,  end=49933052),
  list(name="USP7",     chr="chr16", start=8840635,   end=8940934),
  list(name="KDM1A",    chr="chr1",  start=23035569,  end=23135868),
  list(name="PAX5",     chr="chr9",  start=36992067,  end=37092366),
  list(name="CAMK2A",   chr="chr5",  start=150164623, end=150264922),
  list(name="EPHB1",    chr="chr3",  start=135041420, end=135141719),
  list(name="ZDHHC22",  chr="chr14", start=77092523,  end=77192822)
)

cts <- c("ASO_VPA_vs_Scramble_VPA","ASO_VPA_vs_ASO_CTRL",
         "ASO_CTRL_vs_Scramble_CTRL","Scramble_VPA_vs_Scramble_CTRL")

cat("gene\tchr\tzoom_start\tzoom_end\tbest_dmr_start\tbest_dmr_end\tpval\tcontrast\n")

for (g in genes) {
  gr <- GRanges(g$chr, IRanges(g$start, g$end))
  best_p <- 1
  best_info <- NULL
  for (ct in cts) {
    rds <- paste0("results/dmr/dmr_", ct, ".rds")
    d <- subsetByOverlaps(readRDS(rds), gr)
    if (length(d) > 0) {
      b <- d[which.min(d$pValue)]
      if (b$pValue < best_p) {
        best_p <- b$pValue
        ctr <- as.integer((start(b)+end(b))/2)
        best_info <- paste(c(g$name, g$chr,
          ctr-15000, ctr+15000,
          start(b), end(b),
          formatC(b$pValue, format="e", digits=2),
          ct), collapse="\t")
      }
    }
  }
  if (!is.null(best_info)) cat(best_info, "\n")
  else cat(g$name, "\tNO_DMRs\n")
}
