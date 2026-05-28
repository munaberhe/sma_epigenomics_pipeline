#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --output=logs/star_index.log
#SBATCH --time=03:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

mkdir -p data/reference/star_index

STAR --runMode genomeGenerate \
     --genomeDir data/reference/star_index \
     --genomeFastaFiles data/reference/hg38.fa \
     --sjdbGTFfile data/reference/hg38.ensGene.gtf \
     --runThreadN 8
