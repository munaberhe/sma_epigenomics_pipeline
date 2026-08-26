#!/bin/bash
#SBATCH --job-name=align
#SBATCH --array=1-12%3
#SBATCH --partition=compute
#SBATCH --constraint=ehc
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/data/home/bt25018/sma_epigenomics_pipeline/logs/align/%A_%a.out
#SBATCH --error=/data/home/bt25018/sma_epigenomics_pipeline/logs/align/%A_%a.err

export PATH=/data/home/bt25018/.conda/envs/sma_epigenomics_pipeline/bin:$PATH
set -euo pipefail

SAMPLES=(
  ASO_CTRL_1 ASO_CTRL_2 ASO_CTRL_3
  ASO_VPA_1  ASO_VPA_2  ASO_VPA_3
  Scramble_CTRL_1 Scramble_CTRL_2 Scramble_CTRL_3
  Scramble_VPA_1  Scramble_VPA_2  Scramble_VPA_3
)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
TRIMMED=/gpfs/scratch/bt25018/sma_epigenomics_pipeline/results/trimmed
OUT=/gpfs/scratch/bt25018/sma_epigenomics_pipeline/results/alignments
REF=/gpfs/scratch/bt25018/sma_epigenomics_pipeline/data/reference

mkdir -p ${OUT}/tmp_${SAMPLE}

echo "[$(date)] Aligning: $SAMPLE"
bismark --genome $REF \
  --parallel 4 \
  -1 ${TRIMMED}/${SAMPLE}_1_val_1.fq.gz \
  -2 ${TRIMMED}/${SAMPLE}_2_val_2.fq.gz \
  --output_dir $OUT \
  --temp_dir ${OUT}/tmp_${SAMPLE}

BAM=$(ls ${OUT}/${SAMPLE}*bismark*bt2_pe.bam 2>/dev/null | head -1 || true)
if [ -n "$BAM" ] && [ -f "$BAM" ]; then
  samtools quickcheck "$BAM" || { echo "BAM failed quickcheck: $BAM"; exit 1; }
  echo "[$(date)] Alignment OK - removing trimmed FASTQs for $SAMPLE"
  rm -f ${TRIMMED}/${SAMPLE}_1_val_1.fq.gz \
        ${TRIMMED}/${SAMPLE}_2_val_2.fq.gz
  rm -rf ${OUT}/tmp_${SAMPLE}
else
  echo "[$(date)] ERROR: BAM not found for $SAMPLE"
  exit 1
fi
echo "[$(date)] Done: $SAMPLE"
