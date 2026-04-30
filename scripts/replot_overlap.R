.libPaths("~/R/library")
library(DMRcaller)

contrasts <- list(
  list(name = "ASO_VPA_vs_ASO_CTRL",
       rds  = "results/dmr/ASO_VPA_vs_ASO_CTRL/ASO_VPA_vs_ASO_CTRL_all_chr.rds",
       out  = "results/dmr/ASO_VPA_vs_ASO_CTRL/ASO_VPA_vs_ASO_CTRL_overlap_profile.pdf"),
  list(name = "VPA_vs_Scramble_CTRL",
       rds  = "results/dmr/VPA_vs_Scramble_CTRL/VPA_vs_Scramble_CTRL_all_chr.rds",
       out  = "results/dmr/VPA_vs_Scramble_CTRL/VPA_vs_Scramble_CTRL_overlap_profile.pdf"),
  list(name = "ASO_VPA_vs_Scramble_CTRL",
       rds  = "results/dmr/ASO_VPA_vs_Scramble_CTRL/ASO_VPA_vs_Scramble_CTRL_all_chr.rds",
       out  = "results/dmr/ASO_VPA_vs_Scramble_CTRL/ASO_VPA_vs_Scramble_CTRL_overlap_profile.pdf")
)

# chromosome lengths for hg38
chr_lengths <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
  chr21=46709983,  chr22=50818468,  chrX=156040895)

for (ct in contrasts) {
  message("\nProcessing: ", ct$name)
  if (!file.exists(ct$rds)) { message("  RDS not found, skipping"); next }

  all_dmrs <- readRDS(ct$rds)
  message("  DMRs loaded: ", length(all_dmrs))

  dmr_chroms <- unique(as.character(seqnames(all_dmrs)))
  overlap_list <- GRangesList()

  for (chr in dmr_chroms) {
    if (!chr %in% names(chr_lengths)) next
    chr_dmrs   <- all_dmrs[seqnames(all_dmrs) == chr]
    if (length(chr_dmrs) == 0) next
    chr_region <- GRanges(chr, IRanges(1, chr_lengths[chr]))
    tryCatch({
      op <- computeOverlapProfile(chr_dmrs, chr_region,
                                  windowSize = 100000, binary = FALSE)
      overlap_list[[chr]] <- op
      message("  ", chr, ": ", sum(op$score > 0), " non-zero windows")
    }, error = function(e) message("  ", chr, " error: ", e$message))
  }

  message("  Plotting ", length(overlap_list), " chromosomes...")
  pdf(ct$out, width = 14, height = 6)
  par(mar = c(4, 4, 3, 1) + 0.1)
  plotOverlapProfile(overlap_list, title = ct$name)
  dev.off()
  message("  Saved: ", ct$out)
}

message("\nDone")
