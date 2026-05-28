# Snakefile - SMA Epigenomics WGBS Pipeline
# Muna Berhe, bt25018, MSc Bioinformatics, QMUL
# Supervisor: Prof Radu Zabet
#
# Full pipeline from raw FASTQ to DMR annotation and SMN locus plots.
# Covers both unmasked (genome-wide DMR calling) and SMN1-masked
# (SMN2-specific locus analysis) alignments.
#
# To run: snakemake --profile configs/slurm_profile --jobs 12
# Dry run: snakemake -n

configfile: "configs/config.yaml"

SAMPLES = [
    "ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3",
    "ASO_VPA_1",  "ASO_VPA_2",  "ASO_VPA_3",
    "Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
    "Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3",
]

CONDITIONS = ["ASO_CTRL", "ASO_VPA", "Scramble_CTRL", "Scramble_VPA"]
CHROMS     = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

# Three contrasts matching the locked DMRcaller parameters (confirmed 5 May 2026)
CONTRASTS = [
    "ASO_VPA_vs_Scramble_CTRL",
    "ASO_CTRL_vs_Scramble_CTRL",
    "Scramble_VPA_vs_Scramble_CTRL",
]

RSCRIPT = "/share/apps/rocky9/containers/R/4.5.1/bin/Rscript"
R_LIBS  = "/data/home/bt25018/R/library"


rule all:
    # Final outputs expected from the complete pipeline run
    input:
        # Unmasked alignment done markers (12 samples)
        expand("results/alignments/bs/.{sample}.align.done", sample=SAMPLES),
        # Per-chromosome CpG reports for genome-wide DMR calling
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
        # SMN1-masked alignment done markers
        expand("results/alignments_smn1_masked/bs/.{sample}.align.done", sample=SAMPLES),
        # Masked CX reports - one per sample, gzipped
        expand("results/alignments_smn1_masked/cx_report/.{sample}.cx.done", sample=SAMPLES),
        # Chr5 splits for SMN2 locus analysis
        expand("results/alignments_smn1_masked/chr5_cx/{sample}_chr5.CX_report.txt", sample=SAMPLES),
        # Genome-wide DMR results for all three contrasts
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
        # ChIPseeker annotation with GO/KEGG enrichment
        expand("results/dmr_annotation/{contrast}_annotated.csv", contrast=CONTRASTS),
        # QC plots - PCA and permutation null
        "results/dmr_qc/sample_PCA.pdf",
        "results/dmr_qc/permutation_null_distribution.pdf",
        # Top 10 hypomethylated candidates per contrast (composite ranking)
        expand("results/dmr_annotation/{contrast}_top10_hypo_v2.csv", contrast=CONTRASTS),
        # Coverage QC - replicates vs pooled curves
        "results/qc/coverage_4lines/coverage_4lines_per_condition.pdf",
        # SMN locus plots - main figures for thesis
        "results/smn2_locus_final/SMN_locus_masked_all_comparisons.pdf",
        "results/smn2_locus_final/SMN_locus_masked_lowres_smooth.pdf",


rule trim_galore:
    # Quality trimming and adapter removal. Q20 threshold, min 20bp post-trim.
    # FastQC reports generated alongside trimmed reads.
    input:
        r1 = "data/raw/{sample}_1.fastq.gz",
        r2 = "data/raw/{sample}_2.fastq.gz",
    output:
        r1  = "data/processed/{sample}_1_val_1.fq.gz",
        r2  = "data/processed/{sample}_2_val_2.fq.gz",
        qc1 = "results/qc/trim/{sample}_1_fastqc.html",
        qc2 = "results/qc/trim/{sample}_2_fastqc.html",
    log:
        "logs/trim_{sample}.log"
    threads: 4
    resources:
        mem_mb  = 8000,
        runtime = 120,
    shell:
        """
        trim_galore --paired --fastqc --cores {threads} \
            --output_dir data/processed/ \
            --fastqc_args "--outdir results/qc/trim/" \
            {input.r1} {input.r2} 2>&1 | tee {log}
        """


rule bismark_align:
    # Bismark alignment to unmasked hg38 using Bowtie2 mode.
    # Parameters: -N 1 (1 mismatch per seed), -L 20 (seed length),
    # score_min L,0,-0.6 (allows ~18bp mismatches in 150bp reads).
    # High memory required - bismark spawns 8 parallel instances x 2 threads.
    input:
        r1  = "data/processed/{sample}_1_val_1.fq.gz",
        r2  = "data/processed/{sample}_2_val_2.fq.gz",
        idx = "data/reference/Bisulfite_Genome",
    output:
        bam    = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        report = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_PE_report.txt",
        done   = touch("results/alignments/bs/.{sample}.align.done"),
    log:
        "logs/align_{sample}.log"
    threads: 24
    resources:
        mem_mb     = 200000,
        runtime    = 1440,
        constraint = "ehc",
    shell:
        """
        bismark --genome data/reference \
            --parallel 8 -p 2 \
            -N 1 -L 20 --score-min L,0,-0.6 \
            -1 {input.r1} -2 {input.r2} \
            -o results/alignments/bs/ \
            --temp_dir results/alignments/bs/tmp_{wildcards.sample} \
            2>&1 | tee {log}
        """


rule dedup_and_extract:
    # PCR deduplication followed by methylation extraction.
    # --CX extracts all cytosine contexts; CHG/CHH deleted after to save space.
    # --cytosine_report generates genome-wide CX report directly (no separate
    # coverage2cytosine step needed).
    input:
        bam  = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        done = "results/alignments/bs/.{sample}.align.done",
    output:
        cx = "results/alignments/bs/{sample}_CX_report.txt.CpG_report.txt.gz",
    log:
        "logs/dedup_{sample}.log"
    threads: 8
    resources:
        mem_mb  = 64000,
        runtime = 1440,
    shell:
        """
        mkdir -p results/alignments/dedup

        deduplicate_bismark --paired --bam \
            --output_dir results/alignments/dedup/ \
            {input.bam} 2>&1 | tee {log}

        DEDUP=results/alignments/dedup/{wildcards.sample}_1_val_1_bismark_bt2_pe.deduplicated.bam

        bismark_methylation_extractor \
            --paired-end --CX --cytosine_report \
            --parallel 8 \
            --genome_folder data/reference \
            --output results/alignments/bs/ \
            "$DEDUP" 2>&1 | tee -a {log}

        # Clean up CHG/CHH files immediately to manage scratch space
        rm -f results/alignments/bs/CHG_context_{wildcards.sample}*.txt
        rm -f results/alignments/bs/CHH_context_{wildcards.sample}*.txt
        """


rule split_by_chr:
    # Split genome-wide CpG report into per-chromosome files.
    # CpG context only - DMRcaller reads these per chromosome for parallelised
    # DMR calling (72 jobs: 24 chromosomes x 3 contrasts).
    input:
        cx = "results/alignments/bs/{sample}_CX_report.txt.CpG_report.txt.gz",
    output:
        done = touch("results/alignments/bs/by_chr/.{sample}.split.done"),
    log:
        "logs/split_{sample}.log"
    shell:
        """
        mkdir -p results/alignments/bs/by_chr
        zcat {input.cx} | awk -v s={wildcards.sample} \
            '$6=="CG" {{print > "results/alignments/bs/by_chr/"s"_"$1".CpG_report.txt"}}' \
            2>&1 | tee {log}
        gzip -f results/alignments/bs/by_chr/{wildcards.sample}_chr*.CpG_report.txt
        """


rule mask_and_index:
    # Mask SMN1 locus with Ns to prevent paralog read misassignment.
    # SMN1 and SMN2 share ~99% sequence identity so reads can align to either.
    # Coordinates: chr5:70,924,941-70,953,015 (hg38, 1-based).
    # BED format is 0-based so we use 70924940.
    # After masking, rebuild Bismark index for the masked reference.
    input:
        ref = "data/reference/hg38.fa",
    output:
        masked_fa = "data/reference_smn1_masked/hg38.fa",
        done      = touch("data/reference_smn1_masked/.index.done"),
    log:
        "logs/smn1_mask_index.log"
    threads: 8
    resources:
        mem_mb  = 32000,
        runtime = 240,
    shell:
        """
        mkdir -p data/reference_smn1_masked
        printf "chr5\t70924940\t70953015\tSMN1\n" \
            > data/reference_smn1_masked/smn1_mask.bed

        bedtools maskfasta \
            -fi {input.ref} \
            -bed data/reference_smn1_masked/smn1_mask.bed \
            -fo {output.masked_fa}

        samtools faidx {output.masked_fa}

        bismark_genome_preparation \
            --parallel {threads} \
            data/reference_smn1_masked/ 2>&1 | tee {log}
        """


rule bismark_align_masked:
    # Bismark alignment to SMN1-masked hg38. Identical parameters to unmasked
    # run so mapping efficiency is directly comparable (~80% in both cases).
    # With SMN1 masked, reads from that region align unambiguously to SMN2.
    input:
        r1  = "data/processed/{sample}_1_val_1.fq.gz",
        r2  = "data/processed/{sample}_2_val_2.fq.gz",
        idx = "data/reference_smn1_masked/.index.done",
    output:
        bam    = "results/alignments_smn1_masked/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        report = "results/alignments_smn1_masked/bs/{sample}_1_val_1_bismark_bt2_PE_report.txt",
        done   = touch("results/alignments_smn1_masked/bs/.{sample}.align.done"),
    log:
        "logs/smn1_masked_align_{sample}.log"
    threads: 24
    resources:
        mem_mb     = 200000,
        runtime    = 2880,
        constraint = "ehc",
    shell:
        """
        bismark --genome data/reference_smn1_masked \
            --parallel 8 -p 2 \
            -N 1 -L 20 --score-min L,0,-0.6 \
            -1 {input.r1} -2 {input.r2} \
            -o results/alignments_smn1_masked/bs/ \
            --temp_dir results/alignments_smn1_masked/bs/tmp_{wildcards.sample} \
            2>&1 | tee {log}
        """


rule dedup_and_extract_masked:
    # Dedup and extraction for masked alignment. Same approach as unmasked.
    # Gzip runs synchronously (no & background) to avoid the race condition
    # where cleanup deleted the source before gzip finished.
    input:
        bam  = "results/alignments_smn1_masked/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        done = "results/alignments_smn1_masked/bs/.{sample}.align.done",
    output:
        cx   = "results/alignments_smn1_masked/cx_report/{sample}.CX_report.txt.gz",
        done = touch("results/alignments_smn1_masked/cx_report/.{sample}.cx.done"),
    log:
        "logs/smn1_masked_dedup_{sample}.log"
    threads: 8
    resources:
        mem_mb  = 64000,
        runtime = 1440,
    shell:
        """
        mkdir -p results/alignments_smn1_masked/dedup \
                 results/alignments_smn1_masked/methylation

        deduplicate_bismark --paired --bam \
            --output_dir results/alignments_smn1_masked/dedup/ \
            {input.bam} 2>&1 | tee {log}

        DEDUP=results/alignments_smn1_masked/dedup/{wildcards.sample}_1_val_1_bismark_bt2_pe.deduplicated.bam

        bismark_methylation_extractor \
            --paired-end --CX --cytosine_report \
            --parallel 8 \
            --genome_folder data/reference_smn1_masked \
            --output results/alignments_smn1_masked/methylation/ \
            "$DEDUP" 2>&1 | tee -a {log}

        # Gzip directly into cx_report - synchronous to prevent race condition
        CX_SRC=$(ls results/alignments_smn1_masked/methylation/{wildcards.sample}*CX_report.txt 2>/dev/null | head -1)
        gzip -c "$CX_SRC" > {output.cx}

        # Clean up large intermediate files once CX report is saved
        rm -f results/alignments_smn1_masked/methylation/CHG_context_{wildcards.sample}*.txt
        rm -f results/alignments_smn1_masked/methylation/CHH_context_{wildcards.sample}*.txt
        rm -f "$CX_SRC"
        rm -f "$DEDUP"
        """


rule split_chr5_masked:
    # Extract chr5 CpG data from masked CX reports for SMN2 locus analysis.
    # Only chr5 needed - SMN1 is at chr5:70.9Mb, SMN2 at chr5:70.0Mb.
    input:
        cx   = "results/alignments_smn1_masked/cx_report/{sample}.CX_report.txt.gz",
        done = "results/alignments_smn1_masked/cx_report/.{sample}.cx.done",
    output:
        chr5 = "results/alignments_smn1_masked/chr5_cx/{sample}_chr5.CX_report.txt",
    log:
        "logs/split_chr5_{sample}.log"
    shell:
        """
        mkdir -p results/alignments_smn1_masked/chr5_cx
        zcat {input.cx} | awk '$1=="chr5"' > {output.chr5} 2>&1 | tee {log}
        """


rule dmr_calling:
    # Genome-wide DMR calling using DMRcaller, per-chromosome approach.
    # Parameters locked with Radu Zabet on 5 May 2026:
    #   bins 300bp, minDiff>=0.20, p<=0.01, minCpGs>=4,
    #   minGap=300bp (prevents infinite merge on VPA contrasts),
    #   score test, CG only, 20 permutation seeds.
    # Parallelised across 24 chromosomes x 3 contrasts = 72 SLURM jobs.
    input:
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
        "results/dmr/dmr_summary.tsv",
    log:
        "logs/dmr_calling.log"
    threads: 4
    resources:
        mem_mb  = 64000,
        runtime = 360,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/06b_dmrcaller_by_chr.R \
            2>&1 | tee {log}
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/06c_dmrcaller_combine_chr.R \
            2>&1 | tee -a {log}
        """


rule dmr_annotate:
    # ChIPseeker annotation of high-confidence DMRs (cytosinesCount >= 6).
    # GO biological process and KEGG pathway enrichment via clusterProfiler.
    # Benjamini-Hochberg correction, hyper and hypo DMRs analysed separately.
    input:
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
    output:
        expand("results/dmr_annotation/{contrast}_annotated.csv", contrast=CONTRASTS),
        expand("results/dmr_annotation/{contrast}_GO_BP_hypo_dotplot.pdf", contrast=CONTRASTS),
        expand("results/dmr_annotation/{contrast}_top10_hypo_v2.csv", contrast=CONTRASTS),
    log:
        "logs/dmr_annotate.log"
    resources:
        mem_mb  = 32000,
        runtime = 120,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/07_dmr_annotate.R \
            2>&1 | tee {log}
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/09_top10_by_methylation_v2.R \
            2>&1 | tee -a {log}
        """


rule dmr_qc:
    # QC plots for DMR signal validation.
    # PCA on condition methylation profiles, sample correlation heatmap,
    # permutation null QQ plots (lambda values confirm signal above noise),
    # DMR size distribution, CpG island context overlap.
    input:
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
    output:
        "results/dmr_qc/sample_PCA.pdf",
        "results/dmr_qc/sample_correlation_heatmap.pdf",
        "results/dmr_qc/permutation_null_distribution.pdf",
        "results/dmr_qc/DMR_size_distribution.pdf",
        "results/dmr_qc/CpG_island_overlap.pdf",
    log:
        "logs/dmr_qc.log"
    resources:
        mem_mb  = 32000,
        runtime = 60,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/08_dmr_qc_analysis.R \
            2>&1 | tee {log}
        """


rule coverage_qc:
    # Coverage retention curves showing per-replicate and pooled depth.
    # Confirms 53.8% CpGs reach >=10x when replicates are pooled vs ~33% each.
    input:
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        "results/qc/coverage_4lines/coverage_4lines_per_condition.pdf",
        "results/qc/coverage_4lines/coverage_4lines_per_condition_summary.tsv",
    log:
        "logs/coverage_qc.log"
    resources:
        mem_mb  = 32000,
        runtime = 120,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/10_coverage_qc.R \
            2>&1 | tee {log}
        """


rule smn_locus_plots:
    # SMN locus methylation profiles for both masked and unmasked alignments.
    # Uses DMRcaller plotLocalMethylationProfile (exon track, mean lines) and
    # a lowres sliding window overview (500bp bins, ggplot2 smooth lines).
    # Main thesis figures showing VPA demethylation effect at SMN2.
    input:
        expand("results/alignments_smn1_masked/chr5_cx/{sample}_chr5.CX_report.txt", sample=SAMPLES),
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        "results/smn2_locus_final/SMN_locus_masked_all_comparisons.pdf",
        "results/smn2_locus_final/SMN_locus_unmasked_all_comparisons.pdf",
        "results/smn2_locus_final/SMN_locus_masked_lowres_smooth.pdf",
        "results/smn2_locus_final/SMN_locus_unmasked_lowres_smooth.pdf",
        "results/smn2_locus_final/SMN_weighted_mean_masked_v2.tsv",
    log:
        "logs/smn_locus_plots.log"
    resources:
        mem_mb  = 128000,
        runtime = 120,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/05_smn2_locus_final.R \
            2>&1 | tee {log}
        """
