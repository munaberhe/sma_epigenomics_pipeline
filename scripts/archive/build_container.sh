#!/bin/bash
#SBATCH --job-name=build_container
#SBATCH --output=logs/build_container.log
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

cd ~/sma_epigenomics_pipeline
git pull origin master

apptainer build sma_pipeline.sif Singularity.def

echo "Container build complete"
