#!/bin/bash
#SBATCH --job-name=genome_prep
#SBATCH --partition=compute
#SBATCH --constraint=ehc
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/data/home/bt25018/sma_epigenomics_pipeline/logs/genome_prep.out
#SBATCH --error=/data/home/bt25018/sma_epigenomics_pipeline/logs/genome_prep.err

export PATH=/data/home/bt25018/.conda/envs/sma_epigenomics_pipeline/bin:$PATH
set -euo pipefail

REF=/gpfs/scratch/bt25018/sma_epigenomics_pipeline/data/reference

echo "[$(date)] Starting genome preparation"
bismark_genome_preparation --parallel 8 $REF
echo "[$(date)] Genome preparation complete"
