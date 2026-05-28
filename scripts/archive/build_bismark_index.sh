#!/bin/bash
#SBATCH --job-name=bismark_index
#SBATCH --output=logs/bismark_index.log
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

bismark_genome_preparation \
    --parallel 8 \
    data/reference/
