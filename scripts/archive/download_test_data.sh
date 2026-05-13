#!/bin/bash
#SBATCH --job-name=download_test
#SBATCH --output=logs/download_test.log
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

mkdir -p data/raw

# Small human RNA-seq sample (paired-end, ~50MB)
fasterq-dump SRR1039508 \
    --outdir data/raw/ \
    --threads 4 \
    --split-files

gzip data/raw/SRR1039508_1.fastq
gzip data/raw/SRR1039508_2.fastq

# Small human WGBS sample (paired-end, ~50MB)
fasterq-dump SRR949193 \
    --outdir data/raw/ \
    --threads 4 \
    --split-files

gzip data/raw/SRR949193_1.fastq
gzip data/raw/SRR949193_2.fastq

echo "Downloads complete"
