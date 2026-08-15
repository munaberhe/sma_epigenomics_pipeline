#!/bin/bash
#SBATCH --job-name=space_mon
#SBATCH --partition=compute
#SBATCH --constraint=ehc
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=72:00:00
#SBATCH --output=/data/home/bt25018/sma_epigenomics_pipeline/logs/space_monitor.out
#SBATCH --error=/data/home/bt25018/sma_epigenomics_pipeline/logs/space_monitor.err

ALIGN_JOB=13997653
EXTRACT_JOB=13997654
THRESHOLD=85

while true; do
  PCT=$(qmquota -s 2>/dev/null | grep scratch | awk '{gsub(/%/,"",$2); print int($2)}')
  echo "[$(date)] Scratch usage: ${PCT}%"

  if [ "$PCT" -ge "$THRESHOLD" ]; then
    echo "[$(date)] WARNING: ${PCT}% used — pausing align and extract jobs"
    scontrol hold $ALIGN_JOB
    scontrol hold $EXTRACT_JOB
    echo "[$(date)] Jobs paused. Manual intervention needed."
    # Send yourself a reminder by writing a flag file
    touch ~/sma_epigenomics_pipeline/logs/SPACE_WARNING_$(date +%Y%m%d_%H%M)
    break
  fi

  sleep 1800  # check every 30 minutes
done
