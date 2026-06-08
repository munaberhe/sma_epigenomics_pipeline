#!/bin/bash
set -euo pipefail

cd /data/scratch/bt25018/sma_epigenomics_pipeline
mkdir -p logs

echo "Submitting Stage 1: mask + index"
JID1=$(sbatch --parsable scripts/01_mask_and_index.sh)
echo "  Stage 1 jobid: $JID1"

echo "Submitting Stage 2: bismark align (depends on $JID1)"
JID2=$(sbatch --parsable --dependency=afterok:$JID1 scripts/02_bismark_align.sh)
echo "  Stage 2 jobid: $JID2"

echo "Submitting Stage 3: dedup + methylation extract (depends on $JID2)"
JID3=$(sbatch --parsable --dependency=afterok:$JID2 scripts/03_dedup_and_extract.sh)
echo "  Stage 3 jobid: $JID3"

echo "Submitting Stage 4: split CX by chromosome (depends on $JID3)"
JID4=$(sbatch --parsable --dependency=afterok:$JID3 scripts/04_split_by_chr.sh)
echo "  Stage 4 jobid: $JID4"

echo
echo "All jobs submitted. Monitor with:"
echo "  squeue -u bt25018"
echo "  tail -f logs/smn1_align_${JID2}_1.log"
echo "  bash scripts/scratch_watchdog.sh"
