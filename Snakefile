configfile: "configs/config.yaml"

SAMPLES = list(config["samples"].keys())

rule all:
    input:
        expand("results/qc/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("results/qc/{sample}_2_fastqc.html", sample=SAMPLES),
        "results/qc/multiqc_report.html",
        expand("data/processed/{sample}_1_val_1.fq.gz", sample=SAMPLES),
        expand("data/processed/{sample}_2_val_2.fq.gz", sample=SAMPLES),
        expand("results/alignments/bs/{sample}_bismark.deduplicated.bam", sample=SAMPLES),
        expand("results/alignments/bs/{sample}_bismark.bismark.cov.gz", sample=SAMPLES)

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
        runtime=120
    threads: 4
    shell:
        "fastqc {input.r1} {input.r2} --outdir results/qc/ --threads {threads} 2> {log}"

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
        runtime=20
    shell:
        "multiqc results/qc/ -o results/qc/ 2> {log}"

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
        runtime=360
    threads: 4
    shell:
        "trim_galore --quality {config[trimming][quality]} "
        "--length {config[trimming][min_length]} "
        "--cores {threads} "
        "--paired --gzip "
        "-o data/processed/ {input.r1} {input.r2} 2> {log}"

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
        mem_mb=32000,
        runtime=2880
    threads: 8
    shell:
        "bismark --bowtie2 "
        "-N 1 "
        "-L 20 "
        "--score_min L,0,-0.6 "
        "--genome {input.index} "
        "--parallel 4 "
        "-o results/alignments/bs/ "
        "-1 {input.r1} -2 {input.r2} 2> {log} && "
        "mv results/alignments/bs/{wildcards.sample}_1_val_1_bismark_bt2_pe.bam {output.bam}"

rule bismark_deduplicate:
    input:
        bam="results/alignments/bs/{sample}_bismark.bam"
    output:
        bam="results/alignments/bs/{sample}_bismark.deduplicated.bam"
    log:
        "logs/bismark_dedup/{sample}.log"
    resources:
        mem_mb=16000,
        runtime=480
    threads: 4
    shell:
        "deduplicate_bismark "
        "--paired "
        "--bam "
        "--output_dir results/alignments/bs/ "
        "{input.bam} 2> {log}"

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
        runtime=1440
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
