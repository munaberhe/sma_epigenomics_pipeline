#!/bin/bash
# ----------------------------------------------------------------------
# Standalone scratch watchdog. Run on the login node from cron, or
# call manually whenever you want a status snapshot.
#
# Usage:
#     bash scripts/scratch_watchdog.sh
#
# Optional cron (every 30 min while pipeline runs):
#     */30 * * * * /data/scratch/bt25018/sma_epigenomics_pipeline/scripts/scratch_watchdog.sh \
#         >> /data/scratch/bt25018/sma_epigenomics_pipeline/logs/scratch_watchdog.log 2>&1
# ----------------------------------------------------------------------

set -u

PROJECT=/data/scratch/bt25018/sma_epigenomics_pipeline
LOG_DIR=$PROJECT/logs
mkdir -p "$LOG_DIR"

WARN_GB=200
CRIT_GB=50

ts() { date '+%Y-%m-%d %H:%M:%S'; }

FREE_KB=$(df -P /data/scratch/bt25018 | awk 'NR==2 {print $4}')
FREE_GB=$(( FREE_KB / 1024 / 1024 ))

echo "[$(ts)] free=${FREE_GB}G  warn=${WARN_GB}G  crit=${CRIT_GB}G"

echo "[$(ts)] top results subdirs:"
du -sh "$PROJECT"/results/* 2>/dev/null | sort -rh | head -10

if (( FREE_GB < CRIT_GB )); then
    echo "[$(ts)] CRITICAL: < ${CRIT_GB} GB free."
    echo "[$(ts)] Suggest deleting in this order:"
    echo "          1. rm -rf $PROJECT/results/alignments_smn1_masked/methylation/*.txt"
    echo "          2. rm -rf $PROJECT/results/alignments_smn1_masked/dedup/*.bam   # only after Stage 4 done"
    echo "          3. rm -rf $PROJECT/results/dmr_0p15 $PROJECT/results/dmr_0p25   # old benchmarks"
    if command -v squeue >/dev/null 2>&1; then
        echo "[$(ts)] Pending jobs that would consume more scratch:"
        squeue -u "$USER" -t PD -h -o "%i %j" || true
    fi
elif (( FREE_GB < WARN_GB )); then
    echo "[$(ts)] WARNING: < ${WARN_GB} GB free, monitor closely."
else
    echo "[$(ts)] OK"
fi
