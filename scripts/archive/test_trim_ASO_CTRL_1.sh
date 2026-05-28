#!/bin/bash
#SBATCH --job-name=test_trim
#SBATCH --output=logs/test_trim_ASO_CTRL_1.log
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

mkdir -p data/processed/

trim_galore \
  --quality 20 \
  --length 20 \
  --cores 4 \
  --paired \
  --gzip \
  -o data/processed/ \
  data/raw/ASO_CTRL_1_1.fastq.gz \
  data/raw/ASO_CTRL_1_2.fastq.gz

echo "Trim Galore test complete for ASO_CTRL_1"
