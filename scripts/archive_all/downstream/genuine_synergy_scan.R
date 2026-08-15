.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")

# ---------------------------------------------------------------------------
# Genome-wide interaction scan
# For every "combo-only" candidate gene (present in combination contrast,
# absent from both single-agent contrasts), pulls RAW per-CpG methylation
# across all four conditions at that gene's DMR window, computes the
# additive prediction, and measures the deviation. Runs a label-swap
# permutation null (200 permutations) on the deviation statistic itself
# to get a proper p-value, not just a raw effect size.
#
# Designed for speed: loads each chromosome ONCE per condition, then tests
# every candidate window on that chromosome from the loaded object,
# avoiding repeated file reads per gene.
# ---------------------------------------------------------------------------

BY_CHR_UNMASK <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/synergy_scan"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

CONDITIONS <- c("ASO_CTRL", "Scramble_CTRL", "ASO_VPA", "Scramble_VPA")
N_PERM <- 200
MIN_READS <- 4

# ---------------------------------------------------------------------------
# Step 1: build the candidate gene list (same logic as find_genuine_synergy.R)
# ---------------------------------------------------------------------------
message("Building candidate list...")
combo <- read.csv("results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv")
aso   <- read.csv("results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv")
vpa   <- read.csv("results/dmr_annotation/Scramble_VPA_vs_Scramble_CTRL_annotated.csv")

best_per_gene <- function(df) {
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  df <- df[order(df$pValue), ]
  df[!duplicated(df$SYMBOL),
     c("SYMBOL","seqnames","start","end","pValue","annotation","distanceToTSS")]
}
combo_genes <- best_per_gene(combo)
aso_genes   <- unique(best_per_gene(aso)$SYMBOL)
vpa_genes   <- unique(best_per_gene(vpa)$SYMBOL)

candidates <- combo_genes[!(combo_genes$SYMBOL %in% aso_genes) &
                          !(combo_genes$SYMBOL %in% vpa_genes), ]
# Exclude uncharacterized loci, pseudogenes, lncRNAs, distal-intergenic calls
candidates <- candidates[!grepl("^LOC|^MIR|^LINC|-AS[0-9]*$|-DT$", candidates$SYMBOL), ]
candidates <- candidates[!grepl("Distal Intergenic", candidates$annotation), ]
message("  ", nrow(candidates), " candidates after filtering (intragenic/promoter only)")

# ---------------------------------------------------------------------------
# Step 2: group candidates by chromosome for efficient loading
# ---------------------------------------------------------------------------
candidates$chr <- candidates$seqnames
chroms <- unique(candidates$chr)
message("  spanning ", length(chroms), " chromosomes")

read_unmasked_cpg <- function(condition, chr) {
  files <- file.path(BY_CHR_UNMASK,
                     sprintf("%s_%d_%s.CpG_report.txt.gz", condition, 1:3, chr))
  files <- files[file.exists(files)]
  if (length(files) == 0) return(NULL)
  grs <- lapply(files, function(f) {
    d <- read.table(gzfile(f), header=FALSE, sep="\t",
                    col.names=c("chr","pos","strand","countM","countU","context","tri"),
                    colClasses=c("character","integer","character","integer",
                                 "integer","character","character"))
    d <- d[d$context=="CG", ]
    GRanges(seqnames=d$chr, ranges=IRanges(d$pos,d$pos), strand=d$strand,
            readsM=d$countM, readsN=d$countM+d$countU,
            context=d$context, trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
}

get_proportion <- function(pooled, region) {
  hits <- subsetByOverlaps(pooled, region)
  hits <- hits[hits$readsN >= MIN_READS]
  if (length(hits) == 0) return(list(prop=NA, n=0, readsM=integer(0), readsN=integer(0)))
  list(prop = sum(hits$readsM) / sum(hits$readsN),
       n = length(hits), readsM = hits$readsM, readsN = hits$readsN)
}

# ---------------------------------------------------------------------------
# Step 3: per-chromosome scan
# ---------------------------------------------------------------------------
CHECKPOINT_DIR <- file.path(OUT_DIR, "checkpoints")
dir.create(CHECKPOINT_DIR, recursive=TRUE, showWarnings=FALSE)

all_results <- list()

for (ch in chroms) {
  checkpoint_file <- file.path(CHECKPOINT_DIR, paste0(ch, "_done.rds"))
  if (file.exists(checkpoint_file)) {
    message("\n=== ", ch, " === already done, loading checkpoint")
    all_results[[length(all_results)+1]] <- readRDS(checkpoint_file)
    next
  }

  message("\n=== ", ch, " ===")
  chr_candidates <- candidates[candidates$chr == ch, ]
  message("  ", nrow(chr_candidates), " candidates on this chromosome")

  message("  loading methylation for all 4 conditions...")
  pooled_by_cond <- setNames(
    lapply(CONDITIONS, function(cond) read_unmasked_cpg(cond, ch)),
    CONDITIONS
  )
  if (any(sapply(pooled_by_cond, is.null))) {
    message("  missing data for this chromosome, skipping")
    saveRDS(data.frame(), checkpoint_file)
    next
  }

  all_results_this_chr <- list()
  for (i in seq_len(nrow(chr_candidates))) {
    g <- chr_candidates[i, ]
    region <- GRanges(ch, IRanges(g$start, g$end))

    obs <- lapply(CONDITIONS, function(cond) get_proportion(pooled_by_cond[[cond]], region))
    names(obs) <- CONDITIONS

    if (any(sapply(obs, function(x) is.na(x$prop)))) next

    aso_effect <- obs$ASO_CTRL$prop - obs$Scramble_CTRL$prop
    vpa_effect <- obs$Scramble_VPA$prop - obs$Scramble_CTRL$prop
    predicted  <- obs$Scramble_CTRL$prop + aso_effect + vpa_effect
    actual     <- obs$ASO_VPA$prop
    deviation  <- actual - predicted

    # ---- permutation null for the deviation statistic ----
    # pool all reads across all 4 conditions at this window, then
    # randomly relabel which pooled reads belong to which condition,
    # recompute the deviation under the null, repeat N_PERM times
    all_M <- c(obs$ASO_CTRL$readsM, obs$Scramble_CTRL$readsM,
               obs$ASO_VPA$readsM, obs$Scramble_VPA$readsM)
    all_N <- c(obs$ASO_CTRL$readsN, obs$Scramble_CTRL$readsN,
               obs$ASO_VPA$readsN, obs$Scramble_VPA$readsN)
    group_sizes <- c(length(obs$ASO_CTRL$readsM), length(obs$Scramble_CTRL$readsM),
                     length(obs$ASO_VPA$readsM), length(obs$Scramble_VPA$readsM))
    n_total <- length(all_M)

    if (n_total < 8 || min(group_sizes) < 2) next  # not enough data to permute meaningfully

    perm_devs <- numeric(N_PERM)
    for (p in seq_len(N_PERM)) {
      idx <- sample(n_total)
      cuts <- cumsum(group_sizes)
      g1 <- idx[1:cuts[1]]
      g2 <- idx[(cuts[1]+1):cuts[2]]
      g3 <- idx[(cuts[2]+1):cuts[3]]
      g4 <- idx[(cuts[3]+1):cuts[4]]
      p_aso_ctrl <- sum(all_M[g1])/sum(all_N[g1])
      p_scr_ctrl <- sum(all_M[g2])/sum(all_N[g2])
      p_aso_vpa  <- sum(all_M[g3])/sum(all_N[g3])
      p_scr_vpa  <- sum(all_M[g4])/sum(all_N[g4])
      perm_aso <- p_aso_ctrl - p_scr_ctrl
      perm_vpa <- p_scr_vpa - p_scr_ctrl
      perm_pred <- p_scr_ctrl + perm_aso + perm_vpa
      perm_devs[p] <- p_aso_vpa - perm_pred
    }
    perm_pval <- mean(abs(perm_devs) >= abs(deviation))

    all_results_this_chr[[length(all_results_this_chr)+1]] <- data.frame(
      SYMBOL=g$SYMBOL, chr=ch, start=g$start, end=g$end,
      annotation=g$annotation, distanceToTSS=g$distanceToTSS,
      ASO_CTRL=obs$ASO_CTRL$prop, Scramble_CTRL=obs$Scramble_CTRL$prop,
      ASO_VPA=obs$ASO_VPA$prop, Scramble_VPA=obs$Scramble_VPA$prop,
      aso_effect=aso_effect, vpa_effect=vpa_effect,
      predicted_combo=predicted, actual_combo=actual,
      deviation=deviation, perm_pval=perm_pval
    )
  }
  # save this chromosome's results as a checkpoint, then add to the running list
  chr_results <- do.call(rbind, all_results_this_chr)
  if (!is.null(chr_results)) {
    saveRDS(chr_results, checkpoint_file)
    all_results[[length(all_results)+1]] <- chr_results
    message("  checkpoint saved: ", checkpoint_file)
  } else {
    saveRDS(data.frame(), checkpoint_file)  # mark as done even if empty
  }

  rm(pooled_by_cond); gc()
}

results_df <- do.call(rbind, all_results)
results_df <- results_df[order(results_df$perm_pval, -abs(results_df$deviation)), ]

out_csv <- file.path(OUT_DIR, "genuine_synergy_scan_results.csv")
write.csv(results_df, out_csv, row.names=FALSE)
message("\n\nSaved: ", out_csv)
message("Total candidates tested: ", nrow(results_df))
message("Candidates with perm_pval < 0.05: ", sum(results_df$perm_pval < 0.05, na.rm=TRUE))

cat("\n=== Top 20 by permutation p-value ===\n")
print(head(results_df[, c("SYMBOL","chr","deviation","perm_pval","aso_effect","vpa_effect",
                          "predicted_combo","actual_combo")], 20), row.names=FALSE)
