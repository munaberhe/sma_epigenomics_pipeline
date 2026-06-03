#!/bin/bash
#SBATCH --job-name=smn_merge_align
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --array=1-12
#SBATCH --output=logs/smn_merge_align_%a_%j.log
#SBATCH --error=logs/smn_merge_align_%a_%j.err

conda activate sma_epigenomics_pipeline 2>/dev/null || true
cd /data/home/bt25018/sma_epigenomics_pipeline

SAMPLES=(ASO_CTRL_1 ASO_CTRL_2 ASO_CTRL_3
         ASO_VPA_1 ASO_VPA_2 ASO_VPA_3
         Scramble_CTRL_1 Scramble_CTRL_2 Scramble_CTRL_3
         Scramble_VPA_1 Scramble_VPA_2 Scramble_VPA_3)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
REFDIR="/data/home/bt25018/sma_epigenomics_pipeline/data/reference_smn_merged"
TRIMDIR="/data/scratch/bt25018/sma_epigenomics_pipeline/data/processed"
OUTDIR="/data/home/bt25018/sma_epigenomics_pipeline/results/alignments_smn_merged"
mkdir -p "$OUTDIR"

echo "Aligning $SAMPLE to SMN merged reference..."

bismark \
  --genome "$REFDIR" \
  --bowtie2 \
  -N 1 -L 20 \
  --score_min L,0,-0.6 \
  --parallel 8 \
  -1 "$TRIMDIR/${SAMPLE}_1_val_1.fq.gz" \
  -2 "$TRIMDIR/${SAMPLE}_2_val_2.fq.gz" \
  --output_dir "$OUTDIR"

echo "Done: $SAMPLE"
