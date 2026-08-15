#!/bin/bash
# submit_extraction_chained.sh
# Submits extraction tasks one at a time but chains them via CHG/CHH cleanup
# trigger rather than waiting for full completion. This allows coverage2cytosine
# (light step) to overlap with the next sample's extraction (heavy step),
# saving ~8h per pair vs running fully sequentially.
#
# Usage: bash scripts/submit_extraction_chained.sh <start_task> <end_task>
# Example: bash scripts/submit_extraction_chained.sh 2 6

START=${1:-2}
END=${2:-6}

PROJECT_DIR=/data/scratch/bt25018/sma_epigenomics_pipeline

echo "Submitting chained extraction tasks $START to $END"

# Submit first task immediately
FIRST_JOB=$(sbatch --parsable \
  --array=${START}%1 --mem=64G --time=4-00:00:00 \
  --cpus-per-task=8 --partition=compute \
  --job-name=smn1_met \
  --output=${PROJECT_DIR}/logs/smn1_methext_%A_%a.log \
  --error=${PROJECT_DIR}/logs/smn1_methext_%A_%a.err \
  ${PROJECT_DIR}/scripts/03_dedup_and_extract.sh)

echo "Submitted task $START as job $FIRST_JOB"

# Submit watcher for each subsequent task
PREV_TASK=$START
for TASK in $(seq $((START + 1)) $END); do
    PREV_SAMPLE=$(awk "NR==$PREV_TASK" << 'SAMPLES'
ASO_CTRL_1
ASO_CTRL_2
ASO_CTRL_3
ASO_VPA_1
ASO_VPA_2
ASO_VPA_3
Scramble_CTRL_1
Scramble_CTRL_2
Scramble_CTRL_3
Scramble_VPA_1
Scramble_VPA_2
Scramble_VPA_3
SAMPLES
)
    WATCH_FILE="${PROJECT_DIR}/results/alignments_smn1_masked/methylation/.${PREV_SAMPLE}.chgchh.done"

    sbatch \
      --job-name=watch_t${TASK} \
      --output=${PROJECT_DIR}/logs/watch_task${TASK}_%j.log \
      --time=12:00:00 --mem=1G --partition=compute \
      --wrap="
        echo \"[\$(date)] Watching for: ${WATCH_FILE}\"
        while [[ ! -f '${WATCH_FILE}' ]]; do sleep 120; done
        echo \"[\$(date)] Trigger detected — submitting task ${TASK}\"
        sbatch --array=${TASK}%1 --mem=64G --time=4-00:00:00 \
          --cpus-per-task=8 --partition=compute \
          --job-name=smn1_met \
          --output=${PROJECT_DIR}/logs/smn1_methext_%A_%a.log \
          --error=${PROJECT_DIR}/logs/smn1_methext_%A_%a.err \
          ${PROJECT_DIR}/scripts/03_dedup_and_extract.sh
        echo \"[\$(date)] Task ${TASK} submitted.\"
      "

    echo "Submitted watcher for task $TASK (watching for $PREV_SAMPLE CHG/CHH done)"
    PREV_TASK=$TASK
done
