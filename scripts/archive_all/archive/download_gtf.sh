#!/bin/bash
#SBATCH --job-name=download_gtf
#SBATCH --output=logs/download_gtf.log
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2

wget -P data/reference/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz
gunzip data/reference/hg38.ensGene.gtf.gz
