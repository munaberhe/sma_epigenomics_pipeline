#!/bin/bash
#SBATCH --job-name=extract
#SBATCH --array=1-12%2
#SBATCH --partition=compute
#SBATCH --constraint=ehc
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/data/home/bt25018/sma_epigenomics_pipeline/logs/extract_predup/%A_%a.out
#SBATCH --error=/data/home/bt25018/sma_epigenomics_pipeline/logs/extract_predup/%A_%a.err

export PATH=/data/home/bt25018/.conda/envs/sma_epigenomics_pipeline/bin:$PATH
set -euo pipefail

SAMPLES=(
  ASO_CTRL_1 ASO_CTRL_2 ASO_CTRL_3
  ASO_VPA_1  ASO_VPA_2  ASO_VPA_3
  Scramble_CTRL_1 Scramble_CTRL_2 Scramble_CTRL_3
  Scramble_VPA_1  Scramble_VPA_2  Scramble_VPA_3
)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
BAMS=/gpfs/scratch/bt25018/sma_epigenomics_pipeline/results/alignments
OUT=/gpfs/scratch/bt25018/sma_epigenomics_pipeline/results/methylation_predup
REF=/gpfs/scratch/bt25018/sma_epigenomics_pipeline/data/reference
DONE_FILE=${OUT}/.${SAMPLE}.extract.done

mkdir -p "$OUT"

# ── CHECKPOINT: skip if already successfully extracted ─────────────────
if [ -f "$DONE_FILE" ]; then
  echo "[$(date)] SKIP: $SAMPLE already extracted (found $DONE_FILE)"
  exit 0
fi

# Also skip if CX report already exists (implicit checkpoint)
if compgen -G "${OUT}/${SAMPLE}*CX_report.txt.gz" > /dev/null; then
  echo "[$(date)] SKIP: $SAMPLE CX report already exists"
  touch "$DONE_FILE"
  exit 0
fi

IN_BAM=$(ls ${BAMS}/${SAMPLE}*bismark*bt2_pe.bam 2>/dev/null | head -1 || true)
if [ -z "$IN_BAM" ] || [ ! -f "$IN_BAM" ]; then
  echo "[$(date)] ERROR: no BAM found for $SAMPLE"
  exit 1
fi

samtools quickcheck "$IN_BAM" || { echo "[$(date)] BAM failed quickcheck"; exit 1; }

echo "[$(date)] Extracting (no dedup): $IN_BAM"
bismark_methylation_extractor \
  --paired-end --comprehensive --CX --cytosine_report \
  --genome_folder "$REF" --parallel 6 --gzip \
  --output "$OUT" "$IN_BAM"

if compgen -G "${OUT}/${SAMPLE}*CX_report.txt.gz" > /dev/null; then
  echo "[$(date)] Extraction OK - removing BAM for $SAMPLE"
  rm -f "$IN_BAM"
  touch "$DONE_FILE"
else
  echo "[$(date)] ERROR: CX report not found - keeping BAM"
  exit 1
fi
echo "[$(date)] Done: $SAMPLE"
