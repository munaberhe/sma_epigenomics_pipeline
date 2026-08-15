combo <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv')
aso   <- read.csv('results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv')
vpa   <- read.csv('results/dmr_annotation/Scramble_VPA_vs_Scramble_CTRL_annotated.csv')
combo_v_vpa <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_VPA_annotated.csv')
combo_v_aso <- read.csv('results/dmr_annotation/ASO_VPA_vs_ASO_CTRL_annotated.csv')

best_per_gene <- function(df) {
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  df <- df[order(df$pValue), ]
  df[!duplicated(df$SYMBOL), c("SYMBOL","pValue","proportion1","proportion2")]
}

combo_genes      <- best_per_gene(combo)
aso_genes        <- unique(best_per_gene(aso)$SYMBOL)
vpa_genes        <- unique(best_per_gene(vpa)$SYMBOL)
combo_v_vpa_genes <- unique(best_per_gene(combo_v_vpa)$SYMBOL)
combo_v_aso_genes <- unique(best_per_gene(combo_v_aso)$SYMBOL)

# Genuine synergy candidate definition:
#   - strong DMR in combination vs Scramble_CTRL
#   - NOT in ASO-alone gene list (or only very weak)
#   - NOT in VPA-alone gene list (or only very weak)
#   - DOES show up when combo is compared against EITHER single agent
#     (combo_v_vpa OR combo_v_aso) -- meaning adding the second drug changes something
synergy_candidates <- combo_genes[
  !(combo_genes$SYMBOL %in% aso_genes) &
  !(combo_genes$SYMBOL %in% vpa_genes) &
  (combo_genes$SYMBOL %in% combo_v_vpa_genes | combo_genes$SYMBOL %in% combo_v_aso_genes),
]
synergy_candidates <- synergy_candidates[!grepl("^LOC|^MIR|^LINC|-AS[0-9]*$|-DT$",
                                                 synergy_candidates$SYMBOL), ]
synergy_candidates <- synergy_candidates[order(synergy_candidates$pValue), ]

cat("=== Genuine synergy candidates ===\n")
cat("(strong in combination, absent from BOTH single agents, AND changes when\n")
cat(" combo is compared directly against at least one single agent)\n\n")
print(head(synergy_candidates, 25), row.names=FALSE)
cat("\nTotal candidates:", nrow(synergy_candidates), "\n")
