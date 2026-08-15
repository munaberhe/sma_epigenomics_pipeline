#!/bin/bash
#SBATCH --job-name=download_hg38
#SBATCH --output=logs/download_hg38.log
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

mkdir -p data/reference

# Download hg38 FASTA
wget -P data/reference/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Download GTF annotation (Ensembl)
wget -P data/reference/ https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

# Decompress both
gunzip data/reference/hg38.fa.gz
gunzip data/reference/Homo_sapiens.GRCh38.109.gtf.gz
