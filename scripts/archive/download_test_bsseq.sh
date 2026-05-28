#!/bin/bash
#SBATCH --job-name=download_bsseq
#SBATCH --output=logs/download_bsseq.log
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

mkdir -p data/raw

# Download small WGBS test files directly from ENA
wget -O data/raw/SRR949193_1.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949193/SRR949193_1.fastq.gz

wget -O data/raw/SRR949193_2.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949193/SRR949193_2.fastq.gz

echo "WGBS download complete"
ls -lh data/raw/
