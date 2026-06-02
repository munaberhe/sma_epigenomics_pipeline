#!/bin/bash
#SBATCH --job-name=smn1_mask_index
#SBATCH --output=logs/smn1_mask_index_%j.log
#SBATCH --error=logs/smn1_mask_index_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=compute
#SBATCH --requeue

# Step 1: mask SMN1 in GRCh38 and build Bismark index on the masked reference.
# SMN1 and SMN2 share >99% sequence identity so Bismark discards ambiguously
# mapping reads by default. Masking SMN1 with Ns forces all SMN-derived reads
# to map to SMN2, giving us full coverage at the SMN2 locus.
# Note: this means we can't distinguish SMN1 vs SMN2 methylation — that's
# acceptable because the question is locus-level, not paralog-resolved.

set -euo pipefail

PROJECT_DIR=/data/scratch/bt25018/sma_epigenomics_pipeline
cd "$PROJECT_DIR"
mkdir -p logs data/reference_smn1_masked

source scripts/check_scratch.sh 50

REF_SRC=data/reference/hg38.fa
REF_DST_DIR=data/reference_smn1_masked
REF_DST=$REF_DST_DIR/hg38.fa
MASK_BED=$REF_DST_DIR/smn1_mask.bed
DONE_INDEX=$REF_DST_DIR/.index.done

source ~/.bashrc
conda activate sma_epigenomics_pipeline

# skip if already done
if [[ -f "$DONE_INDEX" ]]; then
    echo "[$(date)] Index already built, skipping."
    exit 0
fi

# SMN1 coordinates on GRCh38 (0-based BED format)
echo "[$(date)] Writing SMN1 mask BED"
printf "chr5\t70924940\t70953015\tSMN1\n" > "$MASK_BED"

# mask the reference if not already done
if [[ ! -f "$REF_DST" ]]; then
    echo "[$(date)] Masking SMN1: $REF_SRC -> $REF_DST"
    bedtools maskfasta \
        -fi "$REF_SRC" \
        -bed "$MASK_BED" \
        -fo "$REF_DST"
    samtools faidx "$REF_DST"
else
    echo "[$(date)] Masked FASTA already present, reusing."
fi

# quick sanity check — first bases of SMN1 region should all be N
echo "[$(date)] Sanity check (should see Ns):"
samtools faidx "$REF_DST" chr5:70924941-70925000

# build Bismark index on masked reference
echo "[$(date)] Building Bismark index"
bismark_genome_preparation \
    --parallel 8 \
    "$REF_DST_DIR/"

touch "$DONE_INDEX"
echo "[$(date)] Done."
