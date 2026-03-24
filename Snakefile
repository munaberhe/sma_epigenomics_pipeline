configfile: "configs/config.yaml"

SAMPLES = config["samples"]

rule all:
    input:
        # QC
        expand("results/qc/{sample}_fastqc.html", sample=SAMPLES),
        "results/qc/multiqc_report.html",
        # Trimming
        expand("data/processed/{sample}_trimmed.fq.gz", sample=SAMPLES),
        # RNA-seq alignments
        expand("results/alignments/rna/{sample}.bam", sample=SAMPLES),
        # Bismark alignments
        expand("results/alignments/bs/{sample}_bismark.bam", sample=SAMPLES),
        # Methylation extraction
        expand("results/alignments/bs/{sample}_bismark.bismark.cov.gz", sample=SAMPLES),
        # featureCounts
        "results/counts/counts.txt",
        # DESeq2
        "results/differential/deseq2_results.csv",
        "results/figures/volcano_plot.png",
        # DMRcaller
        "results/differential/dmrs.csv",
        # Integration
        "results/differential/dmr_deg_overlap.csv"

rule fastqc:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        html="results/qc/{sample}_fastqc.html",
        zip="results/qc/{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    resources:
        mem_mb=4000,
        runtime=30
    threads: 4
    shell:
        "fastqc {input} --outdir results/qc/ --threads {threads} 2> {log}"

rule multiqc:
    input:
        expand("results/qc/{sample}_fastqc.zip", sample=SAMPLES)
    output:
        "results/qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    resources:
        mem_mb=4000,
        runtime=20
    shell:
        "multiqc results/qc/ -o results/qc/ 2> {log}"

rule trim:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "data/processed/{sample}_trimmed.fq.gz"
    log:
        "logs/trim/{sample}.log"
    resources:
        mem_mb=8000,
        runtime=60
    threads: 4
    shell:
        "trim_galore --quality {config[trimming][quality]} "
        "--length {config[trimming][min_length]} "
        "--cores {threads} "
        "--gzip -o data/processed/ {input} 2> {log}"

rule star_align:
    input:
        reads="data/processed/{sample}_trimmed.fq.gz",
        index=config["star_index"]
    output:
        bam="results/alignments/rna/{sample}.bam"
    log:
        "logs/star/{sample}.log"
    resources:
        mem_mb=40000,
        runtime=120
    threads: 8
    shell:
        "STAR --runThreadN {threads} "
        "--genomeDir {input.index} "
        "--readFilesIn {input.reads} "
        "--readFilesCommand zcat "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix results/alignments/rna/{wildcards.sample}_ "
        "2> {log} && "
        "mv results/alignments/rna/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}"

rule bismark_align:
    input:
        reads="data/processed/{sample}_trimmed.fq.gz",
        index=config["bismark_index"]
    output:
        bam="results/alignments/bs/{sample}_bismark.bam"
    log:
        "logs/bismark/{sample}.log"
    resources:
        mem_mb=32000,
        runtime=360
    threads: 8
    shell:
        "bismark --genome {input.index} "
        "--parallel {threads} "
        "-o results/alignments/bs/ "
        "{input.reads} 2> {log} && "
        "mv results/alignments/bs/{wildcards.sample}_trimmed_bismark_bt2.bam {output.bam}"

rule bismark_extract:
    input:
        bam="results/alignments/bs/{sample}_bismark.bam",
        index=config["bismark_index"]
    output:
        cov="results/alignments/bs/{sample}_bismark.bismark.cov.gz"
    log:
        "logs/bismark_extract/{sample}.log"
    resources:
        mem_mb=16000,
        runtime=240
    threads: 8
    shell:
        "bismark_methylation_extractor "
        "--paired-end "
        "--comprehensive "
        "--CX "
        "--cytosine_report "
        "--genome_folder {input.index} "
        "--parallel {threads} "
        "--gzip "
        "-o results/alignments/bs/ "
        "{input.bam} 2> {log}"

rule featurecounts:
    input:
        bams=expand("results/alignments/rna/{sample}.bam", sample=SAMPLES),
        gtf=config["gtf"]
    output:
        counts="results/counts/counts.txt"
    log:
        "logs/featurecounts.log"
    resources:
        mem_mb=8000,
        runtime=60
    threads: 8
    shell:
        "featureCounts "
        "-T {threads} "
        "-a {input.gtf} "
        "-o {output.counts} "
        "{input.bams} 2> {log}"

rule deseq2:
    input:
        counts="results/counts/counts.txt"
    output:
        results="results/differential/deseq2_results.csv",
        volcano="results/figures/volcano_plot.png"
    log:
        "logs/deseq2.log"
    resources:
        mem_mb=16000,
        runtime=60
    threads: 4
    shell:
        "Rscript scripts/deseq2.R "
        "--counts {input.counts} "
        "--padj {config[deseq2][padj_threshold]} "
        "--lfc {config[deseq2][lfc_threshold]} "
        "--out_results {output.results} "
        "--out_volcano {output.volcano} "
        "2> {log}"

rule dmrcaller:
    input:
        coverage=expand("results/alignments/bs/{sample}_bismark.bismark.cov.gz", sample=SAMPLES)
    output:
        dmrs="results/differential/dmrs.csv"
    log:
        "logs/dmrcaller.log"
    resources:
        mem_mb=16000,
        runtime=120
    threads: 4
    shell:
        "Rscript scripts/dmrcaller.R "
        "--coverage {input.coverage} "
        "--context {config[dmrcaller][context]} "
        "--min_coverage {config[dmrcaller][min_coverage]} "
        "--out {output.dmrs} "
        "2> {log}"

rule integrate:
    input:
        dmrs="results/differential/dmrs.csv",
        degs="results/differential/deseq2_results.csv"
    output:
        overlap="results/differential/dmr_deg_overlap.csv"
    log:
        "logs/integrate.log"
    resources:
        mem_mb=8000,
        runtime=30
    threads: 2
    shell:
        "Rscript scripts/integrate.R "
        "--dmrs {input.dmrs} "
        "--degs {input.degs} "
        "--out {output.overlap} "
        "2> {log}"
