#!/bin/bash
#SBATCH --job-name=test_star
#SBATCH --output=logs/test_star_ASO_CTRL_1.log
#SBATCH --time=04:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

mkdir -p results/alignments/rna

STAR --runThreadN 8 \
     --genomeDir data/reference/star_index \
     --readFilesIn data/processed/ASO_CTRL_1_1_val_1.fq.gz data/processed/ASO_CTRL_1_2_val_2.fq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix results/alignments/rna/ASO_CTRL_1_ \
     --outSAMattributes NH HI AS NM MD

echo "STAR test alignment complete for ASO_CTRL_1"
