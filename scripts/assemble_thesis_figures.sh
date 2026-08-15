#!/bin/bash
# assemble_thesis_figures.sh
# Assembles final thesis figures using ImageMagick convert
# Panel labels added via annotate

SCRATCH="/gpfs/scratch/bt25018/sma_epigenomics_pipeline"
OUT="$SCRATCH/results/thesis_figures/final"
mkdir -p "$OUT"
mkdir -p /tmp/thesis_assembly

# helper: extract one page from PDF as PNG at 300dpi
# usage: extract_page input.pdf page_num output.png
extract_page() {
  convert -density 300 "$1[$2]" -background white -flatten "$3"
}

# helper: add panel label to top-left of image
# usage: add_label input.png "A" output.png
add_label() {
  convert "$1" \
    -gravity NorthWest \
    -font DejaVu-Sans-Bold \
    -pointsize 72 \
    -fill black \
    -annotate +20+20 "$2" \
    "$3"
}

echo "Assembling Fig 5.1 - QC panel..."
extract_page "$SCRATCH/results/thesis_figures/Fig5.1b_global_methylation_violin.pdf" 0 /tmp/thesis_assembly/5.1A.png
extract_page "$SCRATCH/results/thesis_figures/Fig5.1c_correlation_heatmap.pdf" 0 /tmp/thesis_assembly/5.1B.png
extract_page "$SCRATCH/results/thesis_figures/Fig_DMR_sample_heatmap.pdf" 0 /tmp/thesis_assembly/5.1C.png
extract_page "$SCRATCH/results/thesis_figures/Fig_lowres_chr1_1Mb.pdf" 0 /tmp/thesis_assembly/5.1D.png
add_label /tmp/thesis_assembly/5.1A.png "A" /tmp/thesis_assembly/5.1A_lab.png
add_label /tmp/thesis_assembly/5.1B.png "B" /tmp/thesis_assembly/5.1B_lab.png
add_label /tmp/thesis_assembly/5.1C.png "C" /tmp/thesis_assembly/5.1C_lab.png
add_label /tmp/thesis_assembly/5.1D.png "D" /tmp/thesis_assembly/5.1D_lab.png
convert +append /tmp/thesis_assembly/5.1A_lab.png /tmp/thesis_assembly/5.1B_lab.png /tmp/thesis_assembly/top51.png
convert +append /tmp/thesis_assembly/5.1C_lab.png /tmp/thesis_assembly/5.1D_lab.png /tmp/thesis_assembly/bot51.png
convert -append /tmp/thesis_assembly/top51.png /tmp/thesis_assembly/bot51.png "$OUT/Fig5.1_QC_panel.png"
echo "  Done Fig5.1"

echo "Assembling Fig 5.5 - annotation + metagene..."
extract_page "$SCRATCH/results/thesis_figures/Fig_annotation_enrichment_heatmap.pdf" 0 /tmp/thesis_assembly/5.5A.png
extract_page "$SCRATCH/results/thesis_figures/Fig_metagene_profile.pdf" 0 /tmp/thesis_assembly/5.5B.png
add_label /tmp/thesis_assembly/5.5A.png "A" /tmp/thesis_assembly/5.5A_lab.png
add_label /tmp/thesis_assembly/5.5B.png "B" /tmp/thesis_assembly/5.5B_lab.png
convert +append /tmp/thesis_assembly/5.5A_lab.png /tmp/thesis_assembly/5.5B_lab.png "$OUT/Fig5.5_annotation_metagene.png"
echo "  Done Fig5.5"

echo "Assembling Fig 5.9 - candidates (single relevant contrast each)..."
# page 3 = ASO_in_VPA contrast (0-indexed = page 2)
extract_page "$SCRATCH/results/thesis_figures/locus_candidates/01_synergy/RELL2_4contrasts.pdf" 2 /tmp/thesis_assembly/5.9A.png
extract_page "$SCRATCH/results/thesis_figures/locus_candidates/01_synergy/DDIT4L_4contrasts.pdf" 2 /tmp/thesis_assembly/5.9B.png
extract_page "$SCRATCH/results/thesis_figures/locus_candidates/02_ASO_restricted/TMEM179B_4contrasts.pdf" 2 /tmp/thesis_assembly/5.9C.png
extract_page "$SCRATCH/results/thesis_figures/locus_candidates/01_synergy/IRF8_4contrasts.pdf" 2 /tmp/thesis_assembly/5.9D.png
add_label /tmp/thesis_assembly/5.9A.png "A" /tmp/thesis_assembly/5.9A_lab.png
add_label /tmp/thesis_assembly/5.9B.png "B" /tmp/thesis_assembly/5.9B_lab.png
add_label /tmp/thesis_assembly/5.9C.png "C" /tmp/thesis_assembly/5.9C_lab.png
add_label /tmp/thesis_assembly/5.9D.png "D" /tmp/thesis_assembly/5.9D_lab.png
convert +append /tmp/thesis_assembly/5.9A_lab.png /tmp/thesis_assembly/5.9B_lab.png /tmp/thesis_assembly/top59.png
convert +append /tmp/thesis_assembly/5.9C_lab.png /tmp/thesis_assembly/5.9D_lab.png /tmp/thesis_assembly/bot59.png
convert -append /tmp/thesis_assembly/top59.png /tmp/thesis_assembly/bot59.png "$OUT/Fig5.9_candidates.png"
echo "  Done Fig5.9"

echo "Assembling Fig 5.10 - SMN2..."
extract_page "$SCRATCH/results/smn2_locus_final/SMN_locus_masked_ASO_VPA_vs_ASO_CTRL.pdf" 0 /tmp/thesis_assembly/5.10A.png
extract_page "$SCRATCH/results/thesis_figures/smn2_extended_igv/SMN2_extended_IGV_50kb.pdf" 0 /tmp/thesis_assembly/5.10B.png
add_label /tmp/thesis_assembly/5.10A.png "A" /tmp/thesis_assembly/5.10A_lab.png
add_label /tmp/thesis_assembly/5.10B.png "B" /tmp/thesis_assembly/5.10B_lab.png
convert -append /tmp/thesis_assembly/5.10A_lab.png /tmp/thesis_assembly/5.10B_lab.png "$OUT/Fig5.10_SMN2.png"
echo "  Done Fig5.10"

# Single panel figures - just convert to PNG with label
echo "Converting single-panel figures..."
for item in \
  "Fig5.2_DMR_heatmap:$SCRATCH/results/thesis_figures/Fig_DMR_sample_heatmap.pdf" \
  "Fig5.3_diverging_bar:$SCRATCH/results/figures/genomic_distribution/chr_dmr_diverging_4contrasts.pdf" \
  "Fig5.4_volcano:$SCRATCH/results/figures/volcano_plots/volcano_4contrasts.pdf" \
  "Fig5.6_GO:$SCRATCH/results/figures/gokegg_pairwise/GO_4contrasts_combined.pdf" \
  "Fig5.7_KEGG:$SCRATCH/results/figures/gokegg_pairwise/KEGG_4contrasts_combined.pdf" \
  "Fig5.8_UpSet:$SCRATCH/results/figures/upset/upset_4contrasts.pdf"
do
  name="${item%%:*}"
  path="${item##*:}"
  convert -density 300 "$path[0]" -background white -flatten "$OUT/${name}.png"
  echo "  Done $name"
done

# Appendix F - full four-contrast locus plots
echo "Converting appendix F figures..."
for item in \
  "FigF.1_RELL2:$SCRATCH/results/thesis_figures/locus_candidates/01_synergy/RELL2_4contrasts.pdf" \
  "FigF.2_DDIT4L:$SCRATCH/results/thesis_figures/locus_candidates/01_synergy/DDIT4L_4contrasts.pdf" \
  "FigF.3_TMEM179B:$SCRATCH/results/thesis_figures/locus_candidates/02_ASO_restricted/TMEM179B_4contrasts.pdf" \
  "FigF.4_IRF8:$SCRATCH/results/thesis_figures/locus_candidates/01_synergy/IRF8_4contrasts.pdf"
do
  name="${item%%:*}"
  path="${item##*:}"
  convert -density 300 "$path" -background white -flatten "$OUT/${name}.png"
  echo "  Done $name"
done

echo "All done. Output: $OUT"
ls "$OUT"
