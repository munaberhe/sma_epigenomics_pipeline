#!/bin/bash
#SBATCH --job-name=split_chr5
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=compute
#SBATCH --output=logs/split_chr5_%A_%a.log
#SBATCH --error=logs/split_chr5_%A_%a.err

set -euo pipefail

SAMPLES=(
    "ASO_CTRL_1" "ASO_CTRL_2" "ASO_CTRL_3"
    "ASO_VPA_1"  "ASO_VPA_2"  "ASO_VPA_3"
    "Scramble_CTRL_1" "Scramble_CTRL_2" "Scramble_CTRL_3"
    "Scramble_VPA_1"  "Scramble_VPA_2"  "Scramble_VPA_3"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$IDX]}

CX_DIR="results/alignments_smn1_masked/cx_report"
OUT_DIR="results/alignments_smn1_masked/chr5_cx"
mkdir -p "$OUT_DIR"

INPUT="$CX_DIR/${SAMPLE}.CX_report.txt.gz"
OUTPUT="$OUT_DIR/${SAMPLE}_chr5.CX_report.txt"

echo "[$(date)] Extracting chr5 for $SAMPLE"

zcat "$INPUT" | awk '$1=="chr5"' > "$OUTPUT"

echo "[$(date)] Done — $(wc -l < $OUTPUT) CpGs on chr5 for $SAMPLE"
