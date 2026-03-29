#!/bin/bash
#SBATCH --job-name=test_bismark
#SBATCH --output=logs/test_bismark_ASO_CTRL_1.log
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

mkdir -p results/alignments/bs

bismark --genome data/reference \
        --parallel 4 \
        -1 data/processed/ASO_CTRL_1_1_val_1.fq.gz \
        -2 data/processed/ASO_CTRL_1_2_val_2.fq.gz \
        -o results/alignments/bs/

echo "Bismark test alignment complete for ASO_CTRL_1"
