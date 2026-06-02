#!/bin/bash
#SBATCH --job-name=smn1_align
#SBATCH --output=logs/smn1_align_%A_%a.log
#SBATCH --error=logs/smn1_align_%A_%a.err
#SBATCH --time=4-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --constraint=ehc
#SBATCH --array=1-12%2
#SBATCH --requeue

# Step 2: align all 12 samples to the SMN1-masked GRCh38 reference.
# Runs as a SLURM array (12 samples, 2 at a time to manage memory).
# Input: trimmed FASTQ pairs from data/processed/
# Output: BAM files in results/alignments_smn1_masked/bs/
# Bismark parameters follow the lab standard for paired-end WGBS.

set -euo pipefail

PROJECT_DIR=/data/scratch/bt25018/sma_epigenomics_pipeline
cd "$PROJECT_DIR"
mkdir -p logs results/alignments_smn1_masked/bs

source scripts/check_scratch.sh 150

# all 12 samples — array index picks one per job
SAMPLES=(
    "ASO_CTRL_1"  "ASO_CTRL_2"  "ASO_CTRL_3"
    "ASO_VPA_1"   "ASO_VPA_2"   "ASO_VPA_3"
    "Scramble_CTRL_1" "Scramble_CTRL_2" "Scramble_CTRL_3"
    "Scramble_VPA_1"  "Scramble_VPA_2"  "Scramble_VPA_3"
)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$IDX]}

R1=data/processed/${SAMPLE}_1_val_1.fq.gz
R2=data/processed/${SAMPLE}_2_val_2.fq.gz
OUTDIR=results/alignments_smn1_masked/bs
DONE=$OUTDIR/.${SAMPLE}.align.done

module unload spack/0.23.1 2>/dev/null || true
source ~/.bashrc
conda activate sma_epigenomics_pipeline

# skip if already done
if [[ -f "$DONE" ]]; then
    echo "[$(date)] $SAMPLE already aligned, skipping."
    exit 0
fi

# check input files exist
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: Missing trimmed FASTQ for $SAMPLE"
    echo "  R1=$R1"
    echo "  R2=$R2"
    exit 1
fi

echo "[$(date)] Aligning $SAMPLE on $(hostname)"
echo "[$(date)] CPUs: $SLURM_CPUS_PER_TASK, mem: $SLURM_MEM_PER_NODE MB"

# use a temp dir so partial outputs don't pollute the output directory
TMP_OUT=$(mktemp -d -p "$OUTDIR" .${SAMPLE}.tmp.XXXXXX)
trap 'rm -rf "$TMP_OUT"' EXIT

/data/home/bt25018/.conda/envs/sma_epigenomics_pipeline/bin/bismark \
    --genome data/reference_smn1_masked \
    -q \
    -N 1 \
    -L 20 \
    --score_min L,0,-0.6 \
    -p 2 \
    --parallel 16 \
    -X 500 \
    -1 "$R1" \
    -2 "$R2" \
    -o "$TMP_OUT/" \
    --temp_dir "$TMP_OUT/tmp"

# find the output BAM (Bismark naming varies by version)
BAM_OUT=""
for pat in \
    "$TMP_OUT/${SAMPLE}_1_val_1_bismark_bt2_pe.bam" \
    "$TMP_OUT/${SAMPLE}_1_val_1_bismark_bt2_PE.bam"; do
    [[ -f "$pat" ]] && BAM_OUT="$pat" && break
done
[[ -z "$BAM_OUT" ]] && echo "ERROR: BAM not found" && exit 1

mv "$BAM_OUT" "$OUTDIR/"

# move alignment report
for pat in \
    "$TMP_OUT/${SAMPLE}_1_val_1_bismark_bt2_PE_report.txt" \
    "$TMP_OUT/${SAMPLE}_1_val_1_bismark_bt2_pe_report.txt"; do
    [[ -f "$pat" ]] && mv "$pat" "$OUTDIR/" && break
done

touch "$DONE"
echo "[$(date)] $SAMPLE done."
