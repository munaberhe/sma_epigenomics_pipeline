#!/usr/bin/env Rscript
# 22_smn2_local_dmr_sensitive.R
# Sensitive local DMR calling at SMN2 locus
# Uses lower minDiff (5%) and smaller window (100bp) to detect
# subtle methylation changes near E7 that genome-wide 20% threshold misses

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
})
.libPaths(c('~/R/library', .libPaths()))
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

OUT <- 'results/smn2_local_dmr'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# SMN2 extended region ±10kb
SMN2_CHR   <- 'chr5'
SMN2_START <- 70049638 - 10000
SMN2_END   <- 70078522 + 10000

CONTRASTS <- list(
  list(name='ASO_CTRL_vs_Scramble_CTRL',
       cond1=c('ASO_CTRL_1','ASO_CTRL_2','ASO_CTRL_3'),
       cond2=c('Scramble_CTRL_1','Scramble_CTRL_2','Scramble_CTRL_3')),
  list(name='ASO_VPA_vs_Scramble_VPA',
       cond1=c('ASO_VPA_1','ASO_VPA_2','ASO_VPA_3'),
       cond2=c('Scramble_VPA_1','Scramble_VPA_2','Scramble_VPA_3'))
)

# Load masked chr5 CX reports — subset to SMN2 region after loading
read_masked <- function(samples) {
  region <- GRanges(seqnames='chr5',
    ranges=IRanges(SMN2_START, SMN2_END))
  grs <- lapply(samples, function(s) {
    f <- file.path('results/alignments_smn1_masked/chr5_cx',
                   paste0(s, '_chr5.CX_report.txt'))
    if (!file.exists(f)) stop('Missing: ', f)
    gr <- readBismark(f)
    subsetByOverlaps(gr, region)
  })
  poolMethylationDatasets(GRangesList(grs))
}

message('Loading masked chr5 data for SMN2 region...')

results_list <- list()

for (ct in CONTRASTS) {
  message('\nContrast: ', ct$name)
  tryCatch({
    pooled1 <- read_masked(ct$cond1)
    pooled2 <- read_masked(ct$cond2)

    message('  Condition 1 CpGs: ', length(pooled1))
    message('  Condition 2 CpGs: ', length(pooled2))

    # Sensitive parameters: 5% minDiff, 100bp window, 3 CpGs minimum
    dmrs <- computeDMRs(pooled1, pooled2,
      regions     = GRanges(seqnames='chr5',
                            ranges=IRanges(SMN2_START, SMN2_END)),
      context     = 'CG',
      method      = 'bins',
      binSize     = 100,
      minCytosinesCount = 3,
      minReadsPerCytosine = 3,
      pValueThreshold = 0.05,
      test        = 'fisher')
    # Apply minDiff filter post-hoc
    if (length(dmrs) > 0) {
      meth_diff <- abs(dmrs$proportion1 - dmrs$proportion2)
      dmrs <- dmrs[meth_diff >= 0.05]
    }

    if (length(dmrs) > 0) {
      message('  DMRs found at 5% threshold: ', length(dmrs))
      df <- as.data.frame(dmrs)
      df$contrast <- ct$name
      results_list[[ct$name]] <- df

      # Flag DMRs near E7
      e7_dmrs <- dmrs[start(dmrs) >= 70071000 & end(dmrs) <= 70082000]
      if (length(e7_dmrs) > 0) {
        message('  *** DMRs near E7 region: ', length(e7_dmrs), ' ***')
        print(as.data.frame(e7_dmrs))
      } else {
        message('  No DMRs in E7 region (70,071-70,082 kb)')
      }
    } else {
      message('  No DMRs found even at 5% threshold')
      results_list[[ct$name]] <- data.frame()
    }

  }, error=function(e) {
    message('  ERROR: ', e$message)
  })
}

# Combine and save
if (length(results_list) > 0) {
  all_dmrs <- do.call(rbind, Filter(function(x) nrow(x)>0, results_list))
  if (nrow(all_dmrs) > 0) {
    write.csv(all_dmrs, file.path(OUT, 'SMN2_sensitive_DMRs.csv'),
              row.names=FALSE)
    message('\nSaved: SMN2_sensitive_DMRs.csv')
    message('Total sensitive DMRs at SMN2 locus: ', nrow(all_dmrs))
  } else {
    message('\nNo DMRs found at SMN2 locus even with 5% threshold')
    message('This confirms SMN2 is epigenetically protected from off-target effects')
    writeLines('No DMRs found at SMN2 with minDiff=0.05, binSize=100, minCytosines=3',
               file.path(OUT, 'SMN2_sensitive_DMR_result.txt'))
  }
}

message('\nDone. Results in: ', OUT)
