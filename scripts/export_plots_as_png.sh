#!/bin/bash
# export_plots_as_png.sh
# Converts all thesis PDFs to PNG for review
# Originals are never touched -- all PNGs go to results/png_exports/
# Uses ImageMagick convert at 150 DPI (good for screen review)
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics

OUT_DIR="results/png_exports"
mkdir -p "$OUT_DIR"

# Find all PDFs in results/ excluding raw data and RDS files
PDFS=$(find results/ \
  -name "*.pdf" \
  -not -path "*/by_chr/*" \
  -not -path "*/archive/*" \
  | sort)

TOTAL=$(echo "$PDFS" | wc -l)
COUNT=0

echo "Found $TOTAL PDFs to convert..."

for pdf in $PDFS; do
  COUNT=$((COUNT + 1))

  # Build output path mirroring the source structure
  # e.g. results/dmr/plots/foo.pdf -> results/png_exports/dmr/plots/foo.png
  rel_path="${pdf#results/}"
  rel_dir=$(dirname "$rel_path")
  base=$(basename "$pdf" .pdf)

  mkdir -p "$OUT_DIR/$rel_dir"
  out_png="$OUT_DIR/$rel_dir/${base}.png"

  echo "[$COUNT/$TOTAL] $rel_path"

  # -density 150: good for screen review, not too large
  # -flatten: handles transparent backgrounds
  # [0]: take only first page for multi-page PDFs
  convert -density 150 \
          -background white \
          -flatten \
          "${pdf}[0]" \
          "$out_png" 2>/dev/null

done

echo ""
echo "Done. PNGs saved to: $OUT_DIR"
find "$OUT_DIR" -name "*.png" | wc -l
