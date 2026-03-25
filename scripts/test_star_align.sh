#!/bin/bash
#SBATCH --job-name=test_star
#SBATCH --output=logs/test_star.log
#SBATCH --time=02:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

mkdir -p results/alignments/rna

STAR --runThreadN 8 \
     --genomeDir data/reference/star_index \
     --readFilesIn data/raw/SRR1039508_1.fastq.gz data/raw/SRR1039508_2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix results/alignments/rna/test_ \
     --outSAMattributes NH HI AS NM MD

echo "STAR test alignment complete"
