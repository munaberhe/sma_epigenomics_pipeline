.libPaths(c('~/R/library', .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

SMN2_START <- 70049638 - 10000
SMN2_END   <- 70078522 + 10000

read_masked <- function(samples) {
  grs <- lapply(samples, function(s) {
    f <- file.path('results/alignments_smn1_masked/chr5_cx',
                   paste0(s, '_chr5.CX_report.txt'))
    d <- read.table(f, header=FALSE, sep='\t',
      col.names=c('chr','pos','strand','countM','countU','context','tri'),
      colClasses=c('character','integer','character','integer',
                   'integer','character','character'))
    d <- d[d$context=='CG' & d$pos >= SMN2_START & d$pos <= SMN2_END, ]
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos, d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

message("Loading ASO_VPA vs Scramble_VPA...")
p1 <- read_masked(c('ASO_VPA_1','ASO_VPA_2','ASO_VPA_3'))
p2 <- read_masked(c('Scramble_VPA_1','Scramble_VPA_2','Scramble_VPA_3'))

dmrs <- computeDMRs(p1, p2,
  regions=GRanges('chr5', IRanges(SMN2_START, SMN2_END)),
  context='CG', method='bins', binSize=50,
  minCytosinesCount=2, minReadsPerCytosine=3,
  pValueThreshold=0.05, test='score')

meth_diff <- abs(dmrs$proportion1 - dmrs$proportion2)
dmrs <- dmrs[meth_diff >= 0.02]

message("DMRs found: ", length(dmrs))
print(as.data.frame(dmrs))

# Are they near E7 (70076521-70076574)?
e7_dmrs <- dmrs[end(dmrs) >= 70070000 & start(dmrs) <= 70082000]
message("Near E7: ", length(e7_dmrs))
if(length(e7_dmrs) > 0) print(as.data.frame(e7_dmrs))
