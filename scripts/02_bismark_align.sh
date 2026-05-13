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

set -euo pipefail

PROJECT_DIR=/data/scratch/bt25018/sma_epigenomics_pipeline
cd "$PROJECT_DIR"

mkdir -p logs
mkdir -p results/alignments_smn1_masked/bs

source scripts/check_scratch.sh 150

SAMPLES=(
    "ASO_CTRL_1"
    "ASO_CTRL_2"
    "ASO_CTRL_3"
    "ASO_VPA_1"
    "ASO_VPA_2"
    "ASO_VPA_3"
    "Scramble_CTRL_1"
    "Scramble_CTRL_2"
    "Scramble_CTRL_3"
    "Scramble_VPA_1"
    "Scramble_VPA_2"
    "Scramble_VPA_3"
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

if [[ -f "$DONE" ]]; then
    echo "[$(date)] $SAMPLE already aligned, skipping."
    exit 0
fi

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: Missing trimmed FASTQ for $SAMPLE"
    echo "  R1=$R1"
    echo "  R2=$R2"
    exit 1
fi

echo "[$(date)] Aligning $SAMPLE on node $(hostname)"
echo "[$(date)] CPUs allocated: $SLURM_CPUS_PER_TASK, mem: $SLURM_MEM_PER_NODE MB"

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

BAM_OUT=""
for pat in \
    "$TMP_OUT/${SAMPLE}_1_val_1_bismark_bt2_pe.bam" \
    "$TMP_OUT/${SAMPLE}_1_val_1_bismark_bt2_PE.bam"; do
    [[ -f "$pat" ]] && BAM_OUT="$pat" && break
done
[[ -z "$BAM_OUT" ]] && echo "ERROR: BAM not found after alignment" && exit 1
mv "$BAM_OUT" "$OUTDIR/"

for pat in \
    "$TMP_OUT/${SAMPLE}_1_val_1_bismark_bt2_PE_report.txt" \
    "$TMP_OUT/${SAMPLE}_1_val_1_bismark_bt2_pe_report.txt"; do
    [[ -f "$pat" ]] && mv "$pat" "$OUTDIR/" && break
done

touch "$DONE"
echo "[$(date)] $SAMPLE alignment complete."
