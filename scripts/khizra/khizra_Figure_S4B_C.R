library(GenomicRanges)
library(gplots)
library(lattice)

load("/storage/projects/ZabetLab/Khizra/Human_Data_27Sep2022/H9/Processed_datasets/H9_full(width)_200bp_promoters200bp_OL50_03june23.Rda")

H9_results_matrix <- matrix(0, 5, 8)
colnames(H9_results_matrix) <- c("Intergenic", "Promoter", "Intron", "Exon", "5' utr", "3' utr", "intronic TEs", "TEs")
rownames(H9_results_matrix) <- c("Common", "Putative Only", "STARR-seq Only", "Neither", "Whole Genome")

H9_results_matrix[1,1] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "intergenic"))
H9_results_matrix[1,2] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "promoter"))
H9_results_matrix[1,3] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "intron"))
H9_results_matrix[1,4] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "exon"))
H9_results_matrix[1,5] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "five_prime_utr"))
H9_results_matrix[1,6] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "three_prime_utr"))
H9_results_matrix[1,7] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "intronic TEs"))
H9_results_matrix[1,8] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "TEs"))

H9_results_matrix[2,1] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "intergenic"))
H9_results_matrix[2,2] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "promoter"))
H9_results_matrix[2,3] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "intron"))
H9_results_matrix[2,4] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "exon"))
H9_results_matrix[2,5] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "five_prime_utr"))
H9_results_matrix[2,6] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "three_prime_utr"))
H9_results_matrix[2,7] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "intronic TEs"))
H9_results_matrix[2,8] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "TEs"))

H9_results_matrix[3,1] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "intergenic"))
H9_results_matrix[3,2] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "promoter"))
H9_results_matrix[3,3] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "intron"))
H9_results_matrix[3,4] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "exon"))
H9_results_matrix[3,5] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "five_prime_utr"))
H9_results_matrix[3,6] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "three_prime_utr"))
H9_results_matrix[3,7] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "intronic TEs"))
H9_results_matrix[3,8] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "TEs"))

H9_results_matrix[4,1] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "intergenic"))
H9_results_matrix[4,2] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "promoter"))
H9_results_matrix[4,3] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "intron"))
H9_results_matrix[4,4] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "exon"))
H9_results_matrix[4,5] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "five_prime_utr"))
H9_results_matrix[4,6] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "three_prime_utr"))
H9_results_matrix[4,7] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "intronic TEs"))
H9_results_matrix[4,8] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "TEs"))

H9_results_matrix[5,1] <- length(which(H9_full$annotations == "intergenic"))
H9_results_matrix[5,2] <- length(which(H9_full$annotations == "promoter"))
H9_results_matrix[5,3] <- length(which(H9_full$annotations == "intron"))
H9_results_matrix[5,4] <- length(which(H9_full$annotations == "exon"))
H9_results_matrix[5,5] <- length(which(H9_full$annotations == "five_prime_utr"))
H9_results_matrix[5,6] <- length(which(H9_full$annotations == "three_prime_utr"))
H9_results_matrix[5,7] <- length(which(H9_full$annotations == "intronic TEs"))
H9_results_matrix[5,8] <- length(which(H9_full$annotations == "TEs"))

H9_hmap <- H9_results_matrix
H9_results_matrix <- t(H9_results_matrix)
H9_results_matrix_aves <- H9_results_matrix

H9_results_matrix_aves[,1] <- H9_results_matrix[,1]/sum(H9_results_matrix[,1])
H9_results_matrix_aves[,2] <- H9_results_matrix[,2]/sum(H9_results_matrix[,2])
H9_results_matrix_aves[,3] <- H9_results_matrix[,3]/sum(H9_results_matrix[,3])
H9_results_matrix_aves[,4] <- H9_results_matrix[,4]/sum(H9_results_matrix[,4])
H9_results_matrix_aves[,5] <- H9_results_matrix[,5]/sum(H9_results_matrix[,5])

H9_hmap_aves <- H9_hmap

H9_hmap_aves[1,] <- H9_hmap[1,]/sum(H9_hmap[1,])
H9_hmap_aves[2,] <- H9_hmap[2,]/sum(H9_hmap[2,])
H9_hmap_aves[3,] <- H9_hmap[3,]/sum(H9_hmap[3,])
H9_hmap_aves[4,] <- H9_hmap[4,]/sum(H9_hmap[4,])
H9_hmap_aves[5,] <- H9_hmap[5,]/sum(H9_hmap[5,])

H9_hmap_comp <- H9_hmap_aves[1:4,]
H9_hmap_comp[1,] <- log2(H9_hmap_comp[1,]/H9_hmap_aves[5,])
H9_hmap_comp[2,] <- log2(H9_hmap_comp[2,]/H9_hmap_aves[5,])
H9_hmap_comp[3,] <- log2(H9_hmap_comp[3,]/H9_hmap_aves[5,])
H9_hmap_comp[4,] <- log2(H9_hmap_comp[4,]/H9_hmap_aves[5,])



# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cols <- cbbPalette[c(1,4,5,6,8,3,7,2)]
custom_at <- seq(-3, 3, by = (8/60))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(60)))

H9_hmap_comp[H9_hmap_comp < -3] <- -3
H9_hmap_comp[H9_hmap_comp > 3] <- 3


pdf("/storage/projects/ZabetLab/Khizra/Human_Data_27Sep2022/H9/Processed_datasets/Plots/Figure_5.13/H1&H9_Annotations(width).pdf",
width = 14, height = 8, pointsize = 14)
par(mar=c(6,4,5,7), mfrow=c(1,1), cex = 1.2)
par(xpd = NA)
x <- barplot(H9_results_matrix_aves, col=cols,
  ylab="Percentage of Results in Region", main="H9 & STARR-seq\nAnnotation Comparisons",
  yaxt = "none", xaxt = "none")
labs <- colnames(H9_results_matrix_aves)
text(cex=1, x=x, y=-0.04, labs, adj = 1, srt=45)
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(6.25,0.75,rownames(H9_results_matrix_aves), fill=cols, bty = "n")

dev.off()


pdf(file = "/storage/projects/ZabetLab/Khizra/Human_Data_27Sep2022/H9/Processed_datasets/Plots/Figure_5.13/H9_Annotation_Heatmap(width).pdf")
  levelplot(t(H9_hmap_comp),
    at = custom_at,
    main = "Predicted Region\nAnnotation Comparisons H9",
    xlab = "log2 Enrichment Difference",
    ylab = "Classification",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()

