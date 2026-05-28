#!/bin/bash
#SBATCH --job-name=smn1_dedup
#SBATCH --output=logs/smn1_dedup_%A_%a.log
#SBATCH --error=logs/smn1_dedup_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compute
#SBATCH --array=11-12%2

set -euo pipefail
PROJECT_DIR=/data/scratch/bt25018/sma_epigenomics_pipeline
cd "$PROJECT_DIR"

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

DEDUP_BAM="$DEDUP_DIR/${SAMPLE}_1_val_1_bismark_bt2_pe.deduplicated.bam"
if [[ -f "$DEDUP_BAM" ]]; then
    echo "[$(date)] Dedup BAM already exists for $SAMPLE, skipping."
    exit 0
fi

module unload spack/0.23.1 2>/dev/null || true
source ~/.bashrc
conda activate sma_epigenomics_pipeline

BAM=""
for cand in \
    "$ALIGN_DIR/${SAMPLE}_1_val_1_bismark_bt2_pe.bam" \
    "$ALIGN_DIR/${SAMPLE}_1_val_1_bismark_bt2_PE.bam"; do
    [[ -f "$cand" ]] && BAM="$cand" && break
done

if [[ -z "$BAM" ]]; then
    echo "ERROR: No aligned BAM found for $SAMPLE"
    exit 1
fi

echo "[$(date)] Deduplicating $SAMPLE"
deduplicate_bismark --paired --bam --output_dir "$DEDUP_DIR" "$BAM"

# Delete the original bs/ BAM immediately after dedup to free space
echo "[$(date)] Deleting original bs/ BAM for $SAMPLE"
rm -f "$BAM"

echo "[$(date)] Dedup complete for $SAMPLE"
