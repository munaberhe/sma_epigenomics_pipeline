#!/bin/bash
#SBATCH --job-name=smn1_dedup
#SBATCH --output=logs/smn1_dedup_%A_%a.log
#SBATCH --error=logs/smn1_dedup_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compute
#SBATCH --array=1-12%2
#SBATCH --requeue
set -euo pipefail

PROJECT_DIR=/data/scratch/bt25018/sma_epigenomics_pipeline
cd "$PROJECT_DIR"
mkdir -p logs

source scripts/check_scratch.sh 100

SAMPLES=(
    "ASO_CTRL_1" "ASO_CTRL_2" "ASO_CTRL_3"
    "ASO_VPA_1"  "ASO_VPA_2"  "ASO_VPA_3"
    "Scramble_CTRL_1" "Scramble_CTRL_2" "Scramble_CTRL_3"
    "Scramble_VPA_1"  "Scramble_VPA_2"  "Scramble_VPA_3"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$IDX]}

ALIGN_DIR=results/alignments_smn1_masked/bs
DEDUP_DIR=results/alignments_smn1_masked/dedup
mkdir -p "$DEDUP_DIR"

DONE=$DEDUP_DIR/.${SAMPLE}.dedup.done

module unload spack/0.23.1 2>/dev/null || true
source ~/.bashrc
conda activate sma_epigenomics_pipeline

if [[ -f "$DONE" ]]; then
    echo "[$(date)] $SAMPLE already deduped, skipping."
    exit 0
fi

BAM=""
for cand in \
    "$ALIGN_DIR/${SAMPLE}_1_val_1_bismark_bt2_pe.bam" \
    "$ALIGN_DIR/${SAMPLE}_1_val_1_bismark_bt2_PE.bam"; do
    [[ -f "$cand" ]] && BAM="$cand" && break
done

if [[ -z "$BAM" ]]; then
    echo "ERROR: No aligned BAM found for $SAMPLE in $ALIGN_DIR"
    exit 1
fi

echo "[$(date)] Deduplicating $SAMPLE"
deduplicate_bismark --paired --bam --output_dir "$DEDUP_DIR" "$BAM"

touch "$DONE"
echo "[$(date)] $SAMPLE dedup complete."
