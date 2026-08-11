#!/bin/bash
#SBATCH --job-name=trim
#SBATCH --array=1-12%6
#SBATCH --partition=compute
#SBATCH --constraint=ehc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data/home/bt25018/sma_epigenomics_pipeline/logs/trim/%A_%a.out
#SBATCH --error=/data/home/bt25018/sma_epigenomics_pipeline/logs/trim/%A_%a.err
set -euo pipefail
source ~/.bashrc && conda activate sma_epigenomics_pipeline

SAMPLES=(
  ASO_CTRL_1 ASO_CTRL_2 ASO_CTRL_3
  ASO_VPA_1  ASO_VPA_2  ASO_VPA_3
  Scramble_CTRL_1 Scramble_CTRL_2 Scramble_CTRL_3
  Scramble_VPA_1  Scramble_VPA_2  Scramble_VPA_3
)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
RAW=/data/Blizard-ZabetLab/SMA_DNAm
OUT=/data/scratch/bt25018/sma_epigenomics_pipeline/results/trimmed

echo "[$(date)] Trimming: $SAMPLE"
trim_galore --paired --cores 8 \
  ${RAW}/${SAMPLE}_1.fastq.gz \
  ${RAW}/${SAMPLE}_2.fastq.gz \
  --output_dir $OUT
echo "[$(date)] Done: $SAMPLE"
