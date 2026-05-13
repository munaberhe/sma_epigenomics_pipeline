#!/bin/bash
set -euo pipefail

MODE=${1:-dry-run}
PROJECT=/data/scratch/bt25018/sma_epigenomics_pipeline
cd "$PROJECT"

ROOT=results/alignments_smn1_masked

echo "=== Mode: $MODE ==="
echo

SAMPLES=(ASO_CTRL_1 ASO_CTRL_2 ASO_CTRL_3 ASO_VPA_1 ASO_VPA_2 ASO_VPA_3 Scramble_CTRL_1 Scramble_CTRL_2 Scramble_CTRL_3 Scramble_VPA_1 Scramble_VPA_2 Scramble_VPA_3)
MISSING=0
for s in "${SAMPLES[@]}"; do
    if [[ ! -f "$ROOT/by_chr/.${s}.splitchr.done" ]]; then
        echo "MISSING split.done for $s"
        MISSING=$((MISSING + 1))
    fi
done
if (( MISSING > 0 )); then
    echo
    echo "Refusing to clean up — $MISSING samples have not finished Stage 4."
    exit 1
fi
echo "All 12 samples have Stage 4 done markers. OK to proceed."
echo

echo "--- Disk usage before cleanup ---"
du -sh "$ROOT"/* 2>/dev/null
echo

run() {
    if [[ "$MODE" == dry-run ]]; then
        echo "WOULD: $*"
    else
        echo "DOING: $*"
        eval "$@"
    fi
}

echo "--- Per-context methylation files (always safe) ---"
for pat in "CHH_*" "CHG_*" "CpG_OT_*" "CpG_OB_*" "CpG_CTOT_*" "CpG_CTOB_*"; do
    for f in "$ROOT/methylation"/$pat.txt "$ROOT/methylation"/$pat.txt.gz; do
        [[ -f "$f" ]] && run rm -f "$f"
    done
done

echo
echo "--- bedGraph and coverage files ---"
for f in "$ROOT/methylation"/*.bedGraph.gz "$ROOT/methylation"/*.cov.gz; do
    [[ -f "$f" ]] && run rm -f "$f"
done

if [[ "$MODE" == "delete-all" ]]; then
    echo
    echo "--- Dedup BAMs (delete-all mode) ---"
    for f in "$ROOT/dedup"/*.bam; do
        [[ -f "$f" ]] && run rm -f "$f"
    done
    echo
    echo "--- Raw aligned BAMs (delete-all mode) ---"
    for f in "$ROOT/bs"/*.bam; do
        [[ -f "$f" ]] && run rm -f "$f"
    done
else
    echo
    echo "--- Skipping dedup BAMs and raw BAMs (use 'delete-all' to remove) ---"
fi

echo
echo "--- Disk usage after cleanup ---"
du -sh "$ROOT"/* 2>/dev/null
echo
echo "=== Done ($MODE) ==="
