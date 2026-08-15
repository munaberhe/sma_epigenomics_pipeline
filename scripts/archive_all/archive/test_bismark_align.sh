#!/bin/bash
#SBATCH --job-name=test_bismark
#SBATCH --output=logs/test_bismark.log
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

mkdir -p results/alignments/bs

bismark --genome data/reference \
        --parallel 4 \
        -1 data/raw/SRR949193_1.fastq.gz \
        -2 data/raw/SRR949193_2.fastq.gz \
        -o results/alignments/bs/

echo "Bismark test alignment complete"
