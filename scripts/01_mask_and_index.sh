#!/bin/bash
#SBATCH --job-name=smn1_mask_index
#SBATCH --output=logs/smn1_mask_index_%j.log
#SBATCH --error=logs/smn1_mask_index_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=compute
#SBATCH --requeue

set -euo pipefail

PROJECT_DIR=/data/scratch/bt25018/sma_epigenomics_pipeline
cd "$PROJECT_DIR"

mkdir -p logs
mkdir -p data/reference_smn1_masked

source scripts/check_scratch.sh 50

REF_SRC=data/reference/hg38.fa
REF_DST_DIR=data/reference_smn1_masked
REF_DST=$REF_DST_DIR/hg38.fa
MASK_BED=$REF_DST_DIR/smn1_mask.bed
DONE_INDEX=$REF_DST_DIR/.index.done

source ~/.bashrc
conda activate sma_epigenomics_pipeline

if [[ -f "$DONE_INDEX" ]]; then
    echo "[$(date)] Stage 1 already complete, skipping."
    exit 0
fi

echo "[$(date)] Writing SMN1 mask BED"
printf "chr5\t70924940\t70953015\tSMN1\n" > "$MASK_BED"

if [[ ! -f "$REF_DST" ]]; then
    echo "[$(date)] Masking SMN1 in $REF_SRC -> $REF_DST"
    bedtools maskfasta \
        -fi "$REF_SRC" \
        -bed "$MASK_BED" \
        -fo "$REF_DST"
    samtools faidx "$REF_DST"
else
    echo "[$(date)] Masked FASTA already present, reusing."
fi

echo "[$(date)] Sanity check (first 60 bases of SMN1 region should be N):"
samtools faidx "$REF_DST" chr5:70924941-70925000

echo "[$(date)] Running bismark_genome_preparation"
bismark_genome_preparation \
    --parallel 8 \
    "$REF_DST_DIR/"

touch "$DONE_INDEX"
echo "[$(date)] Stage 1 complete."
