#!/bin/bash

METH_DIR="results/alignments_smn1_masked/methylation"
CX_DIR="results/alignments_smn1_masked/cx_report"

declare -A SAMPLES
SAMPLES["logs/smn1_methext_11443324_2.err"]="ASO_CTRL_2"
SAMPLES["logs/smn1_methext_11384514_12.err"]="Scramble_VPA_3"
SAMPLES["logs/smn1_methext_10961105_6.err"]="ASO_VPA_3"
SAMPLES["logs/smn1_methext_11506713_5.err"]="ASO_VPA_2"

declare -A DONE

echo "[$(date)] Rescue watcher started"

while true; do
  for LOG in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$LOG]}"

    if [[ -n "${DONE[$SAMPLE]}" ]]; then
      continue
    fi

    if ls "$CX_DIR/${SAMPLE}.CX_report.txt.gz" 2>/dev/null | grep -q .; then
      echo "[$(date)] $SAMPLE already in cx_report — skipping"
      DONE[$SAMPLE]=1
      continue
    fi

    if grep -q "Finished generating genome-wide cytosine report" "$LOG" 2>/dev/null; then
      echo "[$(date)] DETECTED: $SAMPLE finished — starting rescue gzip"
      SRC=$(ls "$METH_DIR"/${SAMPLE}_1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt 2>/dev/null | head -1)
      if [[ -n "$SRC" ]]; then
        gzip -c "$SRC" > "$CX_DIR/${SAMPLE}.CX_report.txt.gz"
        echo "[$(date)] RESCUED: $SAMPLE -> $CX_DIR/${SAMPLE}.CX_report.txt.gz"
        ls -lh "$CX_DIR/${SAMPLE}.CX_report.txt.gz"
      else
        echo "[$(date)] WARNING: $SAMPLE log shows done but CX_report.txt not found"
      fi
      DONE[$SAMPLE]=1
    fi
  done

  if [[ ${#DONE[@]} -eq 4 ]]; then
    echo "[$(date)] All samples rescued. Watcher exiting."
    break
  fi

  sleep 30
done
