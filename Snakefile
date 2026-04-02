configfile: "configs/config.yaml"

SAMPLES = list(config["samples"].keys())

BASE = "/gpfs/scratch/bt25018/sma_epigenomics_pipeline"

rule all:
    input:
        expand("results/qc/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("results/qc/{sample}_2_fastqc.html", sample=SAMPLES),
        "results/qc/multiqc_report.html",
        expand("data/processed/{sample}_1_val_1.fq.gz", sample=SAMPLES),
        expand("data/processed/{sample}_2_val_2.fq.gz", sample=SAMPLES),
        expand("results/alignments/bs/{sample}_bismark.bam", sample=SAMPLES)

rule fastqc:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"]
    output:
        html1="results/qc/{sample}_1_fastqc.html",
        zip1="results/qc/{sample}_1_fastqc.zip",
        html2="results/qc/{sample}_2_fastqc.html",
        zip2="results/qc/{sample}_2_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    resources:
        mem_mb=4000,
        runtime=360
    threads: 4
    shell:
        "fastqc {BASE}/{input.r1} {BASE}/{input.r2} "
        "--outdir {BASE}/results/qc/ "
        "--threads {threads} "
        "2> {BASE}/{log}"

rule multiqc:
    input:
        expand("results/qc/{sample}_1_fastqc.zip", sample=SAMPLES),
        expand("results/qc/{sample}_2_fastqc.zip", sample=SAMPLES)
    output:
        "results/qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    resources:
        mem_mb=4000,
        runtime=360
    shell:
        "multiqc {BASE}/results/qc/ -o {BASE}/results/qc/ 2> {BASE}/{log}"

rule trim:
    input:
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"]
    output:
        r1="data/processed/{sample}_1_val_1.fq.gz",
        r2="data/processed/{sample}_2_val_2.fq.gz"
    log:
        "logs/trim/{sample}.log"
    resources:
        mem_mb=8000,
        runtime=720
    threads: 4
    shell:
        "cd {BASE}/data/processed && "
        "trim_galore --quality {config[trimming][quality]} "
        "--length {config[trimming][min_length]} "
        "--cores {threads} "
        "--paired --gzip "
        "{BASE}/{input.r1} {BASE}/{input.r2} "
        "2> {BASE}/{log}"

rule bismark_align:
    input:
        r1="data/processed/{sample}_1_val_1.fq.gz",
        r2="data/processed/{sample}_2_val_2.fq.gz",
        index=config["bismark_index"]
    output:
        bam="results/alignments/bs/{sample}_bismark.bam"
    log:
        "logs/bismark/{sample}.log"
    resources:
        mem_mb=128000,
        runtime=4320
    threads: 16
    shell:
        "bismark --bowtie2 "
        "-N 1 "
        "-L 20 "
        "--score_min L,0,-0.6 "
        "--reorder "
        "--mm "
        "--genome {BASE}/{input.index} "
        "--parallel 4 "
        "--temp_dir {BASE}/results/alignments/bs/ "
        "-o {BASE}/results/alignments/bs/ "
        "-1 {BASE}/{input.r1} "
        "-2 {BASE}/{input.r2} "
        "2> {BASE}/{log} && "
        "mv {BASE}/results/alignments/bs/{wildcards.sample}_1_val_1_bismark_bt2_pe.bam "
        "{BASE}/{output.bam}"

rule bismark_deduplicate:
    input:
        bam="results/alignments/bs/{sample}_bismark.bam"
    output:
        bam="results/alignments/bs/{sample}_bismark.deduplicated.bam"
    log:
        "logs/bismark_dedup/{sample}.log"
    resources:
        mem_mb=16000,
        runtime=1440
    threads: 4
    shell:
        "deduplicate_bismark "
        "--paired "
        "--bam "
        "--output_dir {BASE}/results/alignments/bs/ "
        "{BASE}/{input.bam} "
        "2> {BASE}/{log}"

rule bismark_extract:
    input:
        bam="results/alignments/bs/{sample}_bismark.deduplicated.bam",
        index=config["bismark_index"]
    output:
        cov="results/alignments/bs/{sample}_bismark.bismark.cov.gz"
    log:
        "logs/bismark_extract/{sample}.log"
    resources:
        mem_mb=16000,
        runtime=2880
    threads: 8
    shell:
        "bismark_methylation_extractor "
        "--paired-end "
        "--comprehensive "
        "--CX "
        "--cytosine_report "
        "--genome_folder {BASE}/{input.index} "
        "--parallel 4 "
        "--gzip "
        "-o {BASE}/results/alignments/bs/ "
        "{BASE}/{input.bam} "
        "2> {BASE}/{log}"
