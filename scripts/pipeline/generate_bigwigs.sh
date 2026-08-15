#!/bin/bash
#SBATCH --job-name=bigwig
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/data/home/bt25018/sma_epigenomics_pipeline/logs/bigwig.out
#SBATCH --error=/data/home/bt25018/sma_epigenomics_pipeline/logs/bigwig.err

source ~/.bashrc
conda activate sma_epigenomics_pipeline

SCRATCH=/gpfs/scratch/bt25018/sma_epigenomics_pipeline
BY_CHR=$SCRATCH/results/alignments/bs/by_chr
OUT=$SCRATCH/results/bigwigs
mkdir -p $OUT

# chromosome sizes for hg38
CHROM_SIZES=$SCRATCH/data/reference/hg38.chrom.sizes
if [ ! -f $CHROM_SIZES ]; then
  fetchChromSizes hg38 > $CHROM_SIZES
fi

CONDITIONS=(ASO_CTRL Scramble_CTRL ASO_VPA Scramble_VPA)
REPS=(1 2 3)

for COND in "${CONDITIONS[@]}"; do
  echo "Processing $COND..."
  
  # merge all reps and all chromosomes into one bedGraph
  TMPBG=$OUT/${COND}_tmp.bedGraph
  > $TMPBG
  
  for CHR in $(seq 1 22) X Y; do
    CHR="chr${CHR}"
    for REP in "${REPS[@]}"; do
      F=$BY_CHR/${COND}_${REP}_${CHR}.CpG_report.txt.gz
      if [ -f "$F" ]; then
        zcat "$F" | awk '$6=="CG" && ($4+$5)>=3 {
          cov=$4+$5
          meth=$4/cov
          printf "%s\t%d\t%d\t%.4f\n", $1, $2-1, $2, meth
        }' >> $TMPBG
      fi
    done
  done
  
  # sort and convert to bigWig
  sort -k1,1 -k2,2n $TMPBG > $OUT/${COND}_sorted.bedGraph
  bedGraphToBigWig $OUT/${COND}_sorted.bedGraph $CHROM_SIZES $OUT/${COND}_methylation.bw
  rm $TMPBG $OUT/${COND}_sorted.bedGraph
  echo "  Saved: ${COND}_methylation.bw"
done

echo "All done. BigWigs in: $OUT"
