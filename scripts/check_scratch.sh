#!/bin/bash
# ----------------------------------------------------------------------
# Sourced by every stage script as a pre-flight scratch check.
# Usage: source scripts/check_scratch.sh <min_gb_free>
# Exits the calling script with code 1 if scratch is below threshold.
# ----------------------------------------------------------------------

_MIN_GB=${1:-50}

_FREE_KB=$(df -P /data/scratch/bt25018 | awk 'NR==2 {print $4}')
_FREE_GB=$(( _FREE_KB / 1024 / 1024 ))

echo "[$(date)] Scratch free: ${_FREE_GB} GB  (required: ${_MIN_GB} GB)"

if (( _FREE_GB < _MIN_GB )); then
    echo "ERROR: Insufficient scratch space. ${_FREE_GB} GB free, need ${_MIN_GB} GB."
    echo "       Run: bash scripts/scratch_watchdog.sh"
    exit 1
fi
