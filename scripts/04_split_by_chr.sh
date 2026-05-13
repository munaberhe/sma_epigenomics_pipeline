#!/bin/bash
#SBATCH --job-name=smn1_splitchr
#SBATCH --output=logs/smn1_splitchr_%A_%a.log
#SBATCH --error=logs/smn1_splitchr_%A_%a.err
#SBATCH --time=4-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=compute
#SBATCH --array=1-12%6
#SBATCH --requeue

set -euo pipefail

PROJECT_DIR=/data/scratch/bt25018/sma_epigenomics_pipeline
cd "$PROJECT_DIR"

source scripts/check_scratch.sh 30

SAMPLES=(
    "ASO_CTRL_1" "ASO_CTRL_2" "ASO_CTRL_3"
    "ASO_VPA_1"  "ASO_VPA_2"  "ASO_VPA_3"
    "Scramble_CTRL_1" "Scramble_CTRL_2" "Scramble_CTRL_3"
    "Scramble_VPA_1"  "Scramble_VPA_2"  "Scramble_VPA_3"
)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$IDX]}

CX_DIR=results/alignments_smn1_masked/cx_report
OUT_DIR=results/alignments_smn1_masked/by_chr
DONE=$OUT_DIR/.${SAMPLE}.splitchr.done

mkdir -p "$OUT_DIR"

module unload spack/0.23.1 2>/dev/null || true
source ~/.bashrc
conda activate sma_epigenomics_pipeline

if [[ -f "$DONE" ]]; then
    echo "[$(date)] $SAMPLE already split, skipping."
    exit 0
fi

CX_FILE=""
for cand in \
    "$CX_DIR/${SAMPLE}.CX_report.txt.gz" \
    "$CX_DIR/${SAMPLE}.CX_report.txt"; do
    [[ -f "$cand" ]] && CX_FILE="$cand" && break
done
if [[ -z "$CX_FILE" ]]; then
    echo "ERROR: No CX report found for $SAMPLE in $CX_DIR"
    exit 1
fi

echo "[$(date)] Splitting $SAMPLE by chromosome..."

if [[ "$CX_FILE" == *.gz ]]; then
    READ_CMD="zcat $CX_FILE"
else
    READ_CMD="cat $CX_FILE"
fi

$READ_CMD | awk -v sample="$SAMPLE" -v outdir="$OUT_DIR" '
    $6 == "CG" {
        chr = $1
        file = outdir "/" sample "_" chr ".CpG_report.txt"
        print > file
    }
'

for f in "$OUT_DIR/${SAMPLE}_chr"*.CpG_report.txt; do
    [[ -f "$f" ]] && gzip -f "$f"
done

touch "$DONE"
echo "[$(date)] $SAMPLE split complete."
