#!/bin/bash
#SBATCH --job-name=subset_bsseq
#SBATCH --output=logs/subset_bsseq.log
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# Extract first 1 million reads (4 million lines) from each file
zcat data/raw/SRR949193_1.fastq.gz | head -4000000 | gzip > data/raw/test_bsseq_1.fastq.gz
zcat data/raw/SRR949193_2.fastq.gz | head -4000000 | gzip > data/raw/test_bsseq_2.fastq.gz

echo "Subsetting complete"
ls -lh data/raw/test_bsseq*
