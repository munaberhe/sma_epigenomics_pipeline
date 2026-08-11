#!/bin/bash
#SBATCH --job-name=check
#SBATCH --partition=compute
#SBATCH --constraint=ehc
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/data/home/bt25018/sma_epigenomics_pipeline/logs/check_outputs.out
#SBATCH --error=/data/home/bt25018/sma_epigenomics_pipeline/logs/check_outputs.err

SAMPLES=(
  ASO_CTRL_1 ASO_CTRL_2 ASO_CTRL_3
  ASO_VPA_1  ASO_VPA_2  ASO_VPA_3
  Scramble_CTRL_1 Scramble_CTRL_2 Scramble_CTRL_3
  Scramble_VPA_1  Scramble_VPA_2  Scramble_VPA_3
)
OUT=/gpfs/scratch/bt25018/sma_epigenomics_pipeline/results/methylation_predup
PASS=0; FAIL=0

for SAMPLE in "${SAMPLES[@]}"; do
  CX=$(ls ${OUT}/${SAMPLE}*CX_report.txt.gz 2>/dev/null | head -1 || true)
  if [ -n "$CX" ] && [ -f "$CX" ]; then
    SIZE=$(du -sh "$CX" | cut -f1)
    echo "PASS: $SAMPLE — $CX ($SIZE)"
    ((PASS++))
  else
    echo "FAIL: $SAMPLE — no CX report found"
    ((FAIL++))
  fi
done

echo ""
echo "Summary: $PASS/12 passed, $FAIL failed"
