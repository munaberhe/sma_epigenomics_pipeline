.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
dir.create("results/hotspot", showWarnings=FALSE, recursive=TRUE)

chr_mb <- c(chr1=249,chr2=242,chr3=198,chr4=190,chr5=182,chr6=171,
            chr7=159,chr8=145,chr9=138,chr10=134,chr11=135,chr12=133,
            chr13=114,chr14=107,chr15=102,chr16=90,chr17=83,chr18=80,
            chr19=59,chr20=64,chr21=47,chr22=51,chrX=156)
CHRS <- names(chr_mb)

COND_COLOURS <- c("#1F3A5F","#C0392B","#D4A017","#6B7280")

contrasts <- list(
  list(name="ASO_VPA_vs_Scramble_CTRL", rds="results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds"),
  list(name="ASO_VPA_vs_ASO_CTRL",      rds="results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds"),
  list(name="ASO_CTRL_vs_Scramble_CTRL",rds="results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds"),
  list(name="Scramble_VPA_vs_Scramble_CTRL",rds="results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
)

per_chr <- function(ct) {
  gr <- readRDS(ct$rds)
  gr <- gr[gr$context=="CG" & as.character(seqnames(gr)) %in% CHRS]
  df <- data.frame(chr=factor(as.character(seqnames(gr)), levels=CHRS),
                   direction=gr$regionType)
  tbl <- as.data.frame(table(df$chr))
  colnames(tbl) <- c("chr","total")
  tbl$mb <- chr_mb[as.character(tbl$chr)]
  tbl$dmrs_per_mb <- round(tbl$total/tbl$mb, 2)
  tbl$contrast <- ct$name
  tbl
}

all_chr <- do.call(rbind, lapply(contrasts, per_chr))

# Print top 5 per contrast
for (ct in contrasts) {
  sub <- all_chr[all_chr$contrast==ct$name,]
  sub <- sub[order(-sub$dmrs_per_mb),]
  cat("\n===", ct$name, "top 5 by DMRs/Mb ===\n")
  print(head(sub[,c("chr","total","dmrs_per_mb")], 5), row.names=FALSE)
}

# Plot: faceted bar chart per contrast
all_chr$chr <- factor(all_chr$chr, levels=CHRS)
all_chr$contrast_label <- gsub("_vs_", " vs\n", all_chr$contrast)

p <- ggplot(all_chr, aes(x=chr, y=dmrs_per_mb, fill=contrast)) +
  geom_col(width=0.8) +
  scale_fill_manual(values=setNames(COND_COLOURS, sapply(contrasts, `[[`, "name")),
                    name=NULL) +
  facet_wrap(~contrast, ncol=1, scales="free_y") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        strip.text=element_text(face="bold", size=9),
        legend.position="none",
        plot.title=element_text(face="bold")) +
  labs(title="DMR density per chromosome across 4 contrasts",
       subtitle="Identifies consistent autosomal hotspots vs contrast-specific signals",
       x="Chromosome", y="DMRs per Mb")

ggsave("results/hotspot/chr_dmr_density_4contrasts.pdf",
       p, width=12, height=12, device=cairo_pdf)
message("saved: chr_dmr_density_4contrasts.pdf")

write.csv(all_chr, "results/hotspot/chr_dmr_density_all_contrasts.csv",
          row.names=FALSE)
message("Done. Outputs in: results/hotspot/")
