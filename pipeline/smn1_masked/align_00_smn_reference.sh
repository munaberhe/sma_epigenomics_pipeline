#!/bin/bash
#SBATCH --job-name=smn_merge_ref
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --partition=compute
#SBATCH --output=logs/smn_merge_ref_%j.log
#SBATCH --error=logs/smn_merge_ref_%j.err

set -euo pipefail
module load anaconda3 2>/dev/null || true
conda activate sma_epigenomics_pipeline 2>/dev/null || true

REFDIR="data/reference_smn_merged"
mkdir -p "$REFDIR"

echo "Step 1: Extract SMN1 and SMN2 sequences"
samtools faidx data/reference/hg38.fa \
  chr5:70049638-70078522 > "$REFDIR/smn2_sequence.fa"
samtools faidx data/reference/hg38.fa \
  chr5:70924941-70953015 > "$REFDIR/smn1_sequence.fa"

echo "Step 2: Create merged chrSMN with 500N gap between loci"
python3 << 'PYEOF'
import os
ref_dir = "data/reference_smn_merged"

with open(f"{ref_dir}/smn2_sequence.fa") as f:
    lines = f.readlines()
    smn2_seq = "".join(l.strip() for l in lines if not l.startswith(">"))

  
with open(f"{ref_dir}/smn1_sequence.fa") as f:
    lines = f.readlines()
    smn1_seq = "".join(l.strip() for l in lines if not l.startswith(">"))

gap = "N" * 500
merged = smn2_seq + gap + smn1_seq

with open(f"{ref_dir}/chrSMN_merged.fa", "w") as f:
    f.write(">chrSMN_merged\n")
    for i in range(0, len(merged), 60):
        f.write(merged[i:i+60] + "\n")

print(f"SMN2 length: {len(smn2_seq)} bp")
print(f"SMN1 length: {len(smn1_seq)} bp")
print(f"Gap: 500 N")
print(f"Merged chrSMN_merged length: {len(merged)} bp")
PYEOF

echo "Step 3: Create masked hg38 (both SMN1 and SMN2 masked)"
echo -e "chr5\t70049637\t70078522\nchr5\t70924940\t70953015"   > "$REFDIR/smn_mask_regions.bed"

bedtools maskfasta   -fi data/reference/hg38.fa   -bed "$REFDIR/smn_mask_regions.bed"   -fo "$REFDIR/hg38_smn_masked.fa"

samtools faidx "$REFDIR/hg38_smn_masked.fa"
echo "Masked reference created"

echo "Step 4: Append chrSMN_merged to masked reference"
cat "$REFDIR/chrSMN_merged.fa" >> "$REFDIR/hg38_smn_masked.fa"
samtools faidx "$REFDIR/hg38_smn_masked.fa"
echo "Reference ready: $REFDIR/hg38_smn_masked.fa"

echo "Step 5: Build Bismark index"
bismark_genome_preparation \
  --bowtie2 \
  --parallel 4 \
  "$REFDIR/"

echo "Done"
echo "Merged reference index built in: $REFDIR/"
#!/bin/bash
#SBATCH --job-name=smn_index
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=compute
#SBATCH --output=logs/smn_index_%j.log
#SBATCH --error=logs/smn_index_%j.err

set -euo pipefail
conda activate sma_epigenomics_pipeline 2>/dev/null || true

echo "Building Bismark index for SMN merged reference..."
bismark_genome_preparation \
  --bowtie2 \
  --parallel 4 \
  data/reference_smn_merged/

echo "Done"
ls -lh data/reference_smn_merged/Bisulfite_Genome/
