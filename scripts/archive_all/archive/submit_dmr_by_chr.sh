#!/bin/bash
# Submit per-chromosome DMR calling for all 3 contrasts
CHROMS=($(echo chr{1..22} chrX chrY))
CONTRASTS=(
  "ASO_VPA_vs_Scramble_CTRL"
  "Scramble_VPA_vs_Scramble_CTRL"
  "ASO_CTRL_vs_Scramble_CTRL"
)

for contrast in "${CONTRASTS[@]}"; do
  for chr in "${CHROMS[@]}"; do
    sbatch --mem=32G --time=2:00:00 --cpus-per-task=2 \
      --partition=compute \
      --job-name=dmr_${chr} \
      --output=logs/dmr_bychr_${contrast}_${chr}_%j.log \
      --wrap="unset LD_PRELOAD && \
              cd /data/scratch/bt25018/sma_epigenomics_pipeline && \
              R_LIBS_USER=/data/home/bt25018/R/library \
              /share/apps/rocky9/containers/R/4.5.1/bin/Rscript \
              scripts/dmrcaller_by_chr.R ${contrast} ${chr}"
  done
done
echo "Submitted $(( ${#CHROMS[@]} * ${#CONTRASTS[@]} )) jobs"
