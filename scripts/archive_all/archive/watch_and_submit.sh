#!/bin/bash
#SBATCH --job-name=watch_submit
#SBATCH --output=logs/watch_submit_%j.log
#SBATCH --time=12:00:00
#SBATCH --mem=1G
#SBATCH --partition=compute

WATCH_FILE="results/alignments_smn1_masked/methylation/.Scramble_VPA_1.chgchh.done"

echo "[$(date)] Watching for: $WATCH_FILE"

while [[ ! -f "$WATCH_FILE" ]]; do
    sleep 300  # check every 5 minutes
done

echo "[$(date)] CHG/CHH cleanup detected — submitting tasks 2-6"

sbatch --array=2-6%1 --mem=64G --time=4-00:00:00 \
  --cpus-per-task=8 --partition=compute \
  --job-name=smn1_met \
  --output=logs/smn1_methext_%A_%a.log \
  --error=logs/smn1_methext_%A_%a.err \
  scripts/03_dedup_and_extract.sh

echo "[$(date)] Tasks 2-6 submitted."
