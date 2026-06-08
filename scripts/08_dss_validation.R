#!/usr/bin/env Rscript
# 20_dss_replicate_testing.R
# Replicate-level DMR calling using DSS (Dispersion Shrinkage for Sequencing)
# Tests ASO_CTRL vs Scramble_CTRL with n=3 replicates per group
# Complements pooled DMRcaller results — validates publishability

suppressPackageStartupMessages({
  library(DSS)
  library(bsseq)
  library(GenomicRanges)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ChIPseeker)
})
.libPaths(c('~/R/library', .libPaths()))

setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/dss_replicate'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

CHR_DIR <- 'results/alignments/bs/by_chr'
CHRS <- paste0('chr', c(1:22, 'X'))

# Sample groups
ASO_CTRL_samples    <- c('ASO_CTRL_1','ASO_CTRL_2','ASO_CTRL_3')
SCRAMBLE_CTRL_samples <- c('Scramble_CTRL_1','Scramble_CTRL_2','Scramble_CTRL_3')

# Function to load one CpG report for one chromosome
load_cpg <- function(sample, chr) {
  f <- file.path(CHR_DIR, paste0(sample, '_', chr, '.CpG_report.txt.gz'))
  if (!file.exists(f)) return(NULL)
  df <- read.table(gzfile(f), header=FALSE, sep='\t',
    col.names=c('chr','pos','strand','M','U','context','trinuc'))
  df <- df[df$strand=='+',]  # CpG context, forward strand only
  df$cov <- df$M + df$U
  df <- df[df$cov >= 1,]
  data.frame(chr=as.character(df$chr), pos=as.integer(df$pos), N=as.integer(df$cov), X=as.integer(df$M))
}

message('Loading CpG reports for ASO_CTRL vs Scramble_CTRL...')
message('Processing ', length(CHRS), ' chromosomes...')

all_dmrs <- list()

for (chr in CHRS) {
  message('  Processing ', chr, '...')

  # Load all 6 samples for this chromosome
  aso1 <- load_cpg(ASO_CTRL_samples[1], chr)
  aso2 <- load_cpg(ASO_CTRL_samples[2], chr)
  aso3 <- load_cpg(ASO_CTRL_samples[3], chr)
  scr1 <- load_cpg(SCRAMBLE_CTRL_samples[1], chr)
  scr2 <- load_cpg(SCRAMBLE_CTRL_samples[2], chr)
  scr3 <- load_cpg(SCRAMBLE_CTRL_samples[3], chr)

  # Skip if any file missing
  if (any(sapply(list(aso1,aso2,aso3,scr1,scr2,scr3), is.null))) {
    message('    Skipping ', chr, ' — missing files')
    next
  }

  # Skip if too few CpGs
  if (nrow(aso1) < 100) {
    message('    Skipping ', chr, ' — too few CpGs')
    next
  }

  tryCatch({
    # Create BSseq object — merge by position first, then build matrix
    all_pos <- sort(unique(c(aso1$pos, aso2$pos, aso3$pos,
                             scr1$pos, scr2$pos, scr3$pos)))
    chr_vec <- rep(chr, length(all_pos))

    make_vec <- function(d, positions) {
      m <- match(positions, d$pos)
      list(M=ifelse(is.na(m), 0L, d$X[m]),
           N=ifelse(is.na(m), 0L, d$N[m]))
    }

    v1<-make_vec(aso1,all_pos); v2<-make_vec(aso2,all_pos)
    v3<-make_vec(aso3,all_pos); v4<-make_vec(scr1,all_pos)
    v5<-make_vec(scr2,all_pos); v6<-make_vec(scr3,all_pos)

    M_mat <- cbind(v1$M, v2$M, v3$M, v4$M, v5$M, v6$M)
    N_mat <- cbind(v1$N, v2$N, v3$N, v4$N, v5$N, v6$N)
    colnames(M_mat) <- colnames(N_mat) <-
      c(ASO_CTRL_samples, SCRAMBLE_CTRL_samples)

    rm(v1,v2,v3,v4,v5,v6,aso1,aso2,aso3,scr1,scr2,scr3)
    gc()

    bs <- BSseq(chr=chr_vec, pos=all_pos,
                M=M_mat, Cov=N_mat,
                sampleNames=c(ASO_CTRL_samples, SCRAMBLE_CTRL_samples))
    rm(M_mat, N_mat, all_pos, chr_vec); gc()

    # Filter low coverage CpGs (min 3x per sample)
    bs.filtered <- bs[which(rowSums(getCoverage(bs) >= 3) == 6),]

    if (nrow(bs.filtered) < 50) {
      message('    Skipping ', chr, ' — too few covered CpGs after filtering')
      next
    }

    message('    ', chr, ': ', nrow(bs.filtered), ' CpGs at >=3x in all samples')

    # DSS DML test
    dml <- DMLtest(bs.filtered,
      group1 = ASO_CTRL_samples,
      group2 = SCRAMBLE_CTRL_samples,
      smoothing = TRUE,
      smoothing.span = 500)

    # Call DMRs — check DML results first
    sig_dmls <- sum(dml$fdr < 0.05, na.rm=TRUE)
    message("    ", chr, ": ", sig_dmls, " significant DMLs at FDR<0.05")
    if (sig_dmls < 3) {
      message("    Skipping ", chr, " — too few significant DMLs")
      next
    }
    dmrs <- tryCatch(
      callDMR(dml,
        delta       = 0.10,
        p.threshold = 0.05,
        minlen      = 50,
        minCG       = 3,
        dis.merge   = 200),
      error = function(e) {
        message("    callDMR error on ", chr, ": ", e$message)
        return(NULL)
      })

    if (!is.null(dmrs) && nrow(dmrs) > 0) {
      dmrs$chr <- chr
      all_dmrs[[chr]] <- dmrs
      message('    ', chr, ': ', nrow(dmrs), ' DMRs found')
    } else {
      message('    ', chr, ': 0 DMRs')
    }

  }, error = function(e) {
    message('    ERROR on ', chr, ': ', e$message)
  })
}

# Combine all chromosomes
if (length(all_dmrs) == 0) {
  message('No DMRs found — check input files')
  quit(status=1)
}

dmr_df <- do.call(rbind, all_dmrs)
message('\nTotal DSS DMRs: ', nrow(dmr_df))
message('Hypermethylated (diff>0): ', sum(dmr_df$diff.Methy > 0))
message('Hypomethylated (diff<0): ', sum(dmr_df$diff.Methy < 0))

# Save results
write.csv(dmr_df, file.path(OUT, 'DSS_ASO_CTRL_vs_Scramble_CTRL_DMRs.csv'),
  row.names=FALSE)
saveRDS(dmr_df, file.path(OUT, 'DSS_ASO_CTRL_vs_Scramble_CTRL_DMRs.rds'))
message('Saved DSS DMR table')

# Convert to GRanges and annotate
dmr_gr <- GRanges(seqnames=dmr_df$chr,
  ranges=IRanges(start=dmr_df$start, end=dmr_df$end),
  diff=dmr_df$diff.Methy,
  nCG=dmr_df$nCG)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno <- annotatePeak(dmr_gr, tssRegion=c(-2000,2000),
  TxDb=txdb, annoDb='org.Hs.eg.db', verbose=FALSE)
anno_df <- as.data.frame(anno)
write.csv(anno_df, file.path(OUT, 'DSS_ASO_CTRL_annotated.csv'), row.names=FALSE)

# GO enrichment on DSS DMR genes
genes <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
message('Unique genes overlapping DSS DMRs: ', length(genes))

if (length(genes) >= 10) {
  go_res <- enrichGO(gene=genes, OrgDb=org.Hs.eg.db,
    ont='BP', pAdjustMethod='BH', pvalueCutoff=0.05,
    qvalueCutoff=0.2, readable=TRUE)

  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    write.csv(go_res@result, file.path(OUT, 'DSS_GO_BP.csv'), row.names=FALSE)
    message('GO terms found: ', sum(go_res@result$p.adjust < 0.05))

    # Dotplot
    pdf(file.path(OUT, 'DSS_GO_dotplot.pdf'), width=10, height=8)
    print(dotplot(go_res, showCategory=20,
      title='GO BP — DSS replicate-level DMRs\nASO_CTRL vs Scramble_CTRL'))
    dev.off()
    message('Saved GO dotplot')

    # Check for neural terms
    neural_terms <- go_res@result[grep('synap|axon|neuro|cognit|learn|memory',
      go_res@result$Description, ignore.case=TRUE),]
    if (nrow(neural_terms) > 0) {
      message('\n*** NEURAL TERMS IN DSS RESULTS ***')
      print(neural_terms[,c('Description','p.adjust','Count')])
      write.csv(neural_terms, file.path(OUT, 'DSS_neural_GO_terms.csv'),
        row.names=FALSE)
    } else {
      message('No neural GO terms at p.adj<0.05 in DSS results')
    }
  } else {
    message('No significant GO terms')
  }
}

message('\nDone. Results in: ', OUT)
