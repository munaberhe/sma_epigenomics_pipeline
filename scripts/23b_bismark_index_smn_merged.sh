#!/bin/bash
#SBATCH --job-name=smn_index
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=compute
#SBATCH --output=logs/smn_index_%j.log
#SBATCH --error=logs/smn_index_%j.err

set -euo pipefail
conda activate sma_epigenomics_pipeline 2>/dev/null || true

echo "Building Bismark index for SMN merged reference..."
bismark_genome_preparation \
  --bowtie2 \
  --parallel 4 \
  data/reference_smn_merged/

echo "Done"
ls -lh data/reference_smn_merged/Bisulfite_Genome/
