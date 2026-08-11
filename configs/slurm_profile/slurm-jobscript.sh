#!/bin/bash
#SBATCH --job-name={rule}
#SBATCH --output=logs/slurm_{rule}_{jobid}.out
#SBATCH --error=logs/slurm_{rule}_{jobid}.err
#SBATCH --mem={resources.mem_mb}M
#SBATCH --time={resources.runtime}
#SBATCH --cpus-per-task={threads}
#SBATCH --partition=compute

# Load conda environment with pipeline tools
source ~/.bashrc
conda activate sma_epigenomics_pipeline

{exec_job}
