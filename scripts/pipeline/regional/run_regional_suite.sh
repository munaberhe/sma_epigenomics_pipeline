#!/bin/bash
#SBATCH --job-name=regional_suite
#SBATCH --output=/data/home/bt25018/sma_epigenomics_pipeline/scripts/pipeline/regional/logs/regional_suite_%j.out
#SBATCH --error=/data/home/bt25018/sma_epigenomics_pipeline/scripts/pipeline/regional/logs/regional_suite_%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail

export SMA_PROJ="/gpfs/scratch/bt25018/sma_epigenomics_pipeline"
SUITE_DIR="/data/home/bt25018/sma_epigenomics_pipeline/scripts/pipeline/regional"
cd "$SUITE_DIR"

echo "=== 40_regional_hotspot_scan ==="
Rscript 40_regional_hotspot_scan.R

echo "=== 41_smn2_adjacent_interaction ==="
Rscript 41_smn2_adjacent_interaction.R

echo "=== 42_chrX_dmr_composition ==="
Rscript 42_chrX_dmr_composition.R

echo "=== done ==="
ls -lh "${SMA_PROJ}/results/regional_suite"
