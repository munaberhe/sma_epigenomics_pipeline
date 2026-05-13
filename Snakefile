# =============================================================================
# Snakefile — SMA Epigenomics WGBS Pipeline
# Muna Berhe · bt25018 · MSc Bioinformatics, QMUL
# Supervisor: Dr Radu Zabet
#
# Pipeline stages:
#   1. QC + trimming (Trim Galore)
#   2. Original hg38 alignment (Bismark)
#   3. Dedup + methylation extraction (original)
#   4. Split by chromosome (original)
#   5. SMN1 masking + masked index build
#   6. Masked re-alignment (Bismark, SMN1-masked hg38)
#   7. Masked dedup + methylation extraction
#   8. Masked split by chromosome
#   9. DMR calling (DMRcaller, 3 contrasts)
#  10. Coverage QC plots
#  11. SMN locus plots (DMRcaller)
# =============================================================================

configfile: "config.yaml"

SAMPLES = [
    "ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3",
    "ASO_VPA_1",  "ASO_VPA_2",  "ASO_VPA_3",
    "Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
    "Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3",
]

CONDITIONS = ["ASO_CTRL", "ASO_VPA", "Scramble_CTRL", "Scramble_VPA"]

CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

RSCRIPT = "/share/apps/rocky9/containers/R/4.5.1/bin/Rscript"
R_LIBS  = "/data/home/bt25018/R/library"

# ── Top-level targets ────────────────────────────────────────────────────────

rule all:
    input:
        # Original alignment done markers
        expand("results/alignments/bs/.{sample}.align.done", sample=SAMPLES),
        # Coverage QC
        "results/qc/coverage_4lines/coverage_4lines_per_condition.pdf",
        # SMN locus unmasked
        "results/qc/smn_locus/smn_locus_dmrcaller_all_comparisons.pdf",
        # DMR results
        "results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds",
        "results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds",
        "results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds",
        # Masked pipeline
        "data/reference_smn1_masked/.index.done",
        expand("results/alignments_smn1_masked/bs/.{sample}.align.done", sample=SAMPLES),
        expand("results/alignments_smn1_masked/cx_report/.{sample}.cx.done", sample=SAMPLES),
        expand("results/alignments_smn1_masked/by_chr/.{sample}.split.done", sample=SAMPLES),
        # Masked SMN2 plots (produced once masked chr5 reports are ready)
        "results/qc/smn_locus/smn_locus_dmrcaller_all_comparisons_masked.pdf",

# ── Stage 1: Trim Galore QC ──────────────────────────────────────────────────

rule trim_galore:
    input:
        r1 = "data/raw/{sample}_1.fastq.gz",
        r2 = "data/raw/{sample}_2.fastq.gz",
    output:
        r1 = "data/processed/{sample}_1_val_1.fq.gz",
        r2 = "data/processed/{sample}_2_val_2.fq.gz",
        qc1 = "results/qc/trim/{sample}_1_fastqc.html",
        qc2 = "results/qc/trim/{sample}_2_fastqc.html",
    log:
        "logs/trim_{sample}.log"
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120,
    shell:
        """
        trim_galore --paired --fastqc --cores {threads} \
            --output_dir data/processed/ \
            --fastqc_args "--outdir results/qc/trim/" \
            {input.r1} {input.r2} 2>&1 | tee {log}
        """

# ── Stage 2: Original hg38 alignment ─────────────────────────────────────────

rule bismark_align:
    input:
        r1 = "data/processed/{sample}_1_val_1.fq.gz",
        r2 = "data/processed/{sample}_2_val_2.fq.gz",
        idx = "data/reference/Bisulfite_Genome",
    output:
        bam  = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        report = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_PE_report.txt",
        done = touch("results/alignments/bs/.{sample}.align.done"),
    log:
        "logs/align_{sample}.log"
    threads: 24
    resources:
        mem_mb  = 200000,
        runtime = 1440,
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

# ── Stage 3: Dedup + methylation extraction (original) ───────────────────────

rule dedup_and_extract:
    input:
        bam  = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        done = "results/alignments/bs/.{sample}.align.done",
    output:
        dedup_bam = "results/alignments/dedup/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bam",
        cov       = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
    log:
        "logs/dedup_{sample}.log"
    threads: 4
    resources:
        mem_mb  = 24000,
        runtime = 720,
    shell:
        """
        deduplicate_bismark --paired --bam \
            --output_dir results/alignments/dedup/ \
            {input.bam} 2>&1 | tee -a {log}

        bismark_methylation_extractor \
            --paired-end --comprehensive \
            --multicore {threads} \
            --bedGraph --cytosine_report --CX \
            --genome_folder data/reference \
            --output results/alignments/bs/ \
            {output.dedup_bam} 2>&1 | tee -a {log}
        """

# ── Stage 4: Split by chromosome (original) ──────────────────────────────────

rule split_by_chr:
    input:
        cx = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
    output:
        done = touch("results/alignments/bs/by_chr/.{sample}.split.done"),
    log:
        "logs/split_{sample}.log"
    shell:
        """
        mkdir -p results/alignments/bs/by_chr
        zcat {input.cx} | awk -v s={wildcards.sample} \
            '{{print > "results/alignments/bs/by_chr/"s"_"$1".CpG_report.txt"}}' \
            2>&1 | tee {log}
        gzip -f results/alignments/bs/by_chr/{wildcards.sample}_chr*.CpG_report.txt
        """

# ── Stage 5: SMN1 masking + Bismark index ─────────────────────────────────────

rule mask_and_index:
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
        # BED (0-based): chr5:70924940-70953015 = 1-based 70924941-70953015
        printf "chr5\t70924940\t70953015\tSMN1\n" > data/reference_smn1_masked/smn1_mask.bed

        bedtools maskfasta \
            -fi {input.ref} \
            -bed data/reference_smn1_masked/smn1_mask.bed \
            -fo {output.masked_fa}

        samtools faidx {output.masked_fa}

        bismark_genome_preparation \
            --parallel {threads} \
            data/reference_smn1_masked/ 2>&1 | tee {log}
        """

# ── Stage 6: Masked alignment ────────────────────────────────────────────────

rule bismark_align_masked:
    input:
        r1   = "data/processed/{sample}_1_val_1.fq.gz",
        r2   = "data/processed/{sample}_2_val_2.fq.gz",
        idx  = "data/reference_smn1_masked/.index.done",
    output:
        bam    = "results/alignments_smn1_masked/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        report = "results/alignments_smn1_masked/bs/{sample}_1_val_1_bismark_bt2_PE_report.txt",
        done   = touch("results/alignments_smn1_masked/bs/.{sample}.align.done"),
    log:
        "logs/smn1_masked_align_{sample}.log"
    threads: 24
    resources:
        mem_mb  = 200000,
        runtime = 2880,
        constraint = "ehc",
    shell:
        """
        TMP=$(mktemp -d results/alignments_smn1_masked/bs/.{wildcards.sample}.tmp.XXXXXX)
        trap 'rm -rf "$TMP"' EXIT

        bismark --genome data/reference_smn1_masked \
            --parallel 16 -p 2 \
            -N 1 -L 20 --score-min L,0,-0.6 \
            -1 {input.r1} -2 {input.r2} \
            -o "$TMP/" --temp_dir "$TMP/tmp" \
            2>&1 | tee {log}

        mv "$TMP"/{wildcards.sample}_1_val_1_bismark_bt2_pe.bam \
            results/alignments_smn1_masked/bs/ 2>/dev/null || \
        mv "$TMP"/{wildcards.sample}_1_val_1_bismark_bt2_PE.bam \
            results/alignments_smn1_masked/bs/ 2>/dev/null
        mv "$TMP"/{wildcards.sample}_1_val_1_bismark_bt2_PE_report.txt \
            results/alignments_smn1_masked/bs/
        """

# ── Stage 7: Masked dedup + methylation extraction ───────────────────────────

rule dedup_and_extract_masked:
    input:
        bam  = "results/alignments_smn1_masked/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        done = "results/alignments_smn1_masked/bs/.{sample}.align.done",
    output:
        cx_report = "results/alignments_smn1_masked/cx_report/{sample}.CX_report.txt.gz",
        done      = touch("results/alignments_smn1_masked/cx_report/.{sample}.cx.done"),
    log:
        "logs/smn1_masked_dedup_{sample}.log"
    threads: 4
    resources:
        mem_mb  = 24000,
        runtime = 720,
    shell:
        """
        mkdir -p results/alignments_smn1_masked/dedup \
                 results/alignments_smn1_masked/methylation \
                 results/alignments_smn1_masked/cx_report

        deduplicate_bismark --paired --bam \
            --output_dir results/alignments_smn1_masked/dedup/ \
            {input.bam} 2>&1 | tee {log}

        DEDUP=results/alignments_smn1_masked/dedup/{wildcards.sample}_1_val_1_bismark_bt2_pe.deduplicated.bam

        bismark_methylation_extractor \
            --paired-end --comprehensive \
            --multicore {threads} \
            --bedGraph --cytosine_report --CX \
            --genome_folder data/reference_smn1_masked \
            --output results/alignments_smn1_masked/methylation/ \
            "$DEDUP" 2>&1 | tee -a {log}

        CX_SRC=$(ls results/alignments_smn1_masked/methylation/{wildcards.sample}*deduplicated*.CX_report.txt.gz 2>/dev/null | head -1)
        cp -p "$CX_SRC" {output.cx_report}
        """

# ── Stage 8: Masked split by chromosome ──────────────────────────────────────

rule split_by_chr_masked:
    input:
        cx   = "results/alignments_smn1_masked/cx_report/{sample}.CX_report.txt.gz",
        done = "results/alignments_smn1_masked/cx_report/.{sample}.cx.done",
    output:
        done = touch("results/alignments_smn1_masked/by_chr/.{sample}.split.done"),
    log:
        "logs/smn1_masked_split_{sample}.log"
    shell:
        """
        mkdir -p results/alignments_smn1_masked/by_chr
        zcat {input.cx} | awk -v s={wildcards.sample} \
            '{{print > "results/alignments_smn1_masked/by_chr/"s"_"$1".CpG_report.txt"}}' \
            2>&1 | tee {log}
        gzip -f results/alignments_smn1_masked/by_chr/{wildcards.sample}_chr*.CpG_report.txt
        """

# ── Stage 9: DMR calling ──────────────────────────────────────────────────────

rule dmr_calling:
    input:
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        dmr_combo   = "results/dmr/dmr_ASO_VPA_vs_Scramble_CTRL.rds",
        dmr_aso     = "results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds",
        dmr_vpa     = "results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds",
        summary_tsv = "results/dmr/dmr_summary.tsv",
    log:
        "logs/dmr_calling.log"
    threads: 4
    resources:
        mem_mb  = 64000,
        runtime = 360,
    shell:
        """
        module unload miniforge/24.7.1 2>/dev/null || true
        module load R/4.5.1
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/dmrcaller_analysis.R \
            2>&1 | tee {log}
        """

# ── Stage 10: Coverage QC plots ──────────────────────────────────────────────

rule coverage_qc:
    input:
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        pdf = "results/qc/coverage_4lines/coverage_4lines_per_condition.pdf",
        png = "results/qc/coverage_4lines/coverage_4lines_per_condition.png",
        tsv = "results/qc/coverage_4lines/coverage_4lines_per_condition_summary.tsv",
    log:
        "logs/coverage_qc.log"
    resources:
        mem_mb  = 32000,
        runtime = 120,
    shell:
        """
        module unload miniforge/24.7.1 2>/dev/null || true
        module load R/4.5.1
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/coverage_4lines_per_condition.R \
            2>&1 | tee {log}
        """

# ── Stage 11a: SMN locus — unmasked DMRcaller plots ──────────────────────────

rule smn_locus_dmrcaller:
    input:
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        all_pdf  = "results/qc/smn_locus/smn_locus_dmrcaller_all_comparisons.pdf",
        aso_pdf  = "results/qc/smn_locus/smn_locus_dmrcaller_ASO_vs_Scramble_CTRL.pdf",
        vpa_pdf  = "results/qc/smn_locus/smn_locus_dmrcaller_ASO_vs_Scramble_VPA.pdf",
        ctrl_pdf = "results/qc/smn_locus/smn_locus_dmrcaller_VPA_vs_CTRL_ASO.pdf",
        tsv      = "results/qc/smn_locus/smn_locus_dmrcaller_comparisons_summary.tsv",
    log:
        "logs/smn_locus_dmrcaller.log"
    resources:
        mem_mb  = 32000,
        runtime = 120,
    shell:
        """
        module unload miniforge/24.7.1 2>/dev/null || true
        module load R/4.5.1
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/smn_locus_dmrcaller_comparisons.R \
            2>&1 | tee {log}
        """

# ── Stage 11c: SMN locus — masked DMRcaller plots ────────────────────────────

rule smn_locus_dmrcaller_masked:
    input:
        expand("results/alignments_smn1_masked/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        masked_pdf = "results/qc/smn_locus/smn_locus_dmrcaller_all_comparisons_masked.pdf",
    log:
        "logs/smn_locus_dmrcaller_masked.log"
    resources:
        mem_mb  = 32000,
        runtime = 120,
    shell:
        """
        module unload miniforge/24.7.1 2>/dev/null || true
        module load R/4.5.1
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/smn_locus_dmrcaller_comparisons.R \
            --masked TRUE \
            2>&1 | tee {log}
        """
