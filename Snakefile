# Snakefile - SMA Epigenomics WGBS Pipeline
# Muna Berhe, bt25018, MSc Bioinformatics, QMUL
# Supervisor: Prof Radu Zabet

configfile: "configs/config.yaml"

SAMPLES = [
    "ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3",
    "ASO_VPA_1",  "ASO_VPA_2",  "ASO_VPA_3",
    "Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
    "Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3",
]
CONDITIONS = ["ASO_CTRL", "ASO_VPA", "Scramble_CTRL", "Scramble_VPA"]
CHROMS     = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
CONTRASTS  = [
    "ASO_CTRL_vs_Scramble_CTRL",
    "ASO_VPA_vs_Scramble_CTRL",
    "ASO_VPA_vs_ASO_CTRL",
    "ASO_VPA_vs_Scramble_VPA",
    "Scramble_VPA_vs_Scramble_CTRL",
]

RSCRIPT = "/share/apps/rocky9/containers/R/4.5.1/bin/Rscript"
R_LIBS  = "/data/home/bt25018/R/library"

rule all:
    input:
        # Alignment
        expand("results/alignments/bs/.{sample}.align.done", sample=SAMPLES),
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
        expand("results/alignments_smn1_masked/bs/.{sample}.align.done", sample=SAMPLES),
        expand("results/alignments_smn1_masked/chr5_cx/{sample}_chr5.CX_report.txt", sample=SAMPLES),
        # DMR calling
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
        "results/dmr/dmr_summary.tsv",
        # QC
        "results/dmr_qc/sample_PCA_12samples_chr1.pdf",
        "results/qc/coverage_4lines/coverage_4lines_per_condition.pdf",
        "results/qc/additional_qc/mbias_CpG_all_samples.pdf",
        # DMR annotation
        expand("results/dmr_annotation/{contrast}_annotated.csv", contrast=CONTRASTS),
        "results/dmr_annotation/DMR_annotation_combined_count.pdf",
        "results/dmr_overlap/dmr_upset_plot_5contrasts.pdf",
        # Enrichment
        expand("results/dmr_annotation/{contrast}_GO_BP_hypo_dotplot.pdf", contrast=CONTRASTS),
        "results/dmr_annotation/msigdb_v2/ASO_specific_msigdb_all_combined.pdf",
        "results/h3k9me2_overlap/h3k9me2_signal_boxplot.pdf",
        "results/smn2_enhancer/smn2_h3k27ac_peak_summary.csv",
        # TF motif
        "results/tf_motif/motif_enrichment_volcano.pdf",
        # Locus plots
        "results/lowres_profiles/lowres_allgroups_chrX_500kb.pdf",
        "results/lowres_profiles/lowres_allgroups_chr5_SMN2_10kb.pdf",
        expand("results/dmr/plots/annotated/{gene}_annotated.pdf",
               gene=["MTA1-DT","SLC32A1","CHRNB3","GLRA4","GFRA2",
                     "SEMA3C","PHACTR3","SOX5","RNF169","SMN2"]),
        # SMN2 locus
        "results/smn2_locus_final/SMN_locus_masked_all_comparisons.pdf",
        "results/smn2_locus_final/SMN2_masked_vs_unmasked_comparison.pdf",
        "results/smn2_local_dmr/SMN2_sensitive_DMRs_all_contrasts.csv",
        # Validation
        "results/dss_replicate/DSS_GO_dotplot.pdf",
        # Summary plots
        "results/tss_metaplot/TSS_metaplot.pdf",
        "results/tss_metaplot/DMR_heatmap_top500_ASO_VPA_methdiff.pdf",
        "results/manhattan/manhattan_ASO_CTRL_vs_Scramble_CTRL.pdf",
        # Top10
        "results/dmr_annotation/top20_hypo_ASO_by_methdiff.csv",

rule trim_galore:
    input:
        r1 = "data/raw/{sample}_1.fastq.gz",
        r2 = "data/raw/{sample}_2.fastq.gz",
    output:
        r1  = "data/processed/{sample}_1_val_1.fq.gz",
        r2  = "data/processed/{sample}_2_val_2.fq.gz",
        qc1 = "results/qc/trim/{sample}_1_fastqc.html",
        qc2 = "results/qc/trim/{sample}_2_fastqc.html",
    log: "logs/trim_{sample}.log"
    threads: 4
    resources:
        mem_mb=8000, runtime=120,
    shell:
        """
        trim_galore --paired --fastqc --cores {threads} \
            --output_dir data/processed/ \
            --fastqc_args "--outdir results/qc/trim/" \
            {input.r1} {input.r2} 2>&1 | tee {log}
        """

rule bismark_align:
    input:
        r1  = "data/processed/{sample}_1_val_1.fq.gz",
        r2  = "data/processed/{sample}_2_val_2.fq.gz",
        idx = "data/reference/Bisulfite_Genome",
    output:
        bam  = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        done = touch("results/alignments/bs/.{sample}.align.done"),
    log: "logs/align_{sample}.log"
    threads: 24
    resources:
        mem_mb=200000, runtime=1440, constraint="ehc",
    shell:
        """
        bismark --genome data/reference \
            --parallel 8 -p 2 -N 1 -L 20 --score-min L,0,-0.6 \
            -1 {input.r1} -2 {input.r2} \
            -o results/alignments/bs/ \
            --temp_dir results/alignments/bs/tmp_{wildcards.sample} \
            2>&1 | tee {log}
        """

rule dedup_and_extract:
    input:
        bam  = "results/alignments/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        done = "results/alignments/bs/.{sample}.align.done",
    output:
        cx = "results/alignments/bs/{sample}_CX_report.txt.CpG_report.txt.gz",
    log: "logs/dedup_{sample}.log"
    threads: 8
    resources:
        mem_mb=64000, runtime=1440,
    shell:
        """
        mkdir -p results/alignments/dedup
        deduplicate_bismark --paired --bam \
            --output_dir results/alignments/dedup/ {input.bam} 2>&1 | tee {log}
        DEDUP=results/alignments/dedup/{wildcards.sample}_1_val_1_bismark_bt2_pe.deduplicated.bam
        bismark_methylation_extractor \
            --paired-end --CX --cytosine_report --parallel 8 \
            --genome_folder data/reference \
            --output results/alignments/bs/ "$DEDUP" 2>&1 | tee -a {log}
        rm -f results/alignments/bs/CHG_context_{wildcards.sample}*.txt
        rm -f results/alignments/bs/CHH_context_{wildcards.sample}*.txt
        """

rule split_by_chr:
    input:
        cx = "results/alignments/bs/{sample}_CX_report.txt.CpG_report.txt.gz",
    output:
        done = touch("results/alignments/bs/by_chr/.{sample}.split.done"),
    log: "logs/split_{sample}.log"
    shell:
        """
        mkdir -p results/alignments/bs/by_chr
        zcat {input.cx} | awk -v s={wildcards.sample} \
            '$6=="CG" {{print > "results/alignments/bs/by_chr/"s"_"$1".CpG_report.txt"}}' \
            2>&1 | tee {log}
        gzip -f results/alignments/bs/by_chr/{wildcards.sample}_chr*.CpG_report.txt
        """

rule mask_and_index:
    input:
        ref = "data/reference/hg38.fa",
    output:
        masked_fa = "data/reference_smn1_masked/hg38.fa",
        done      = touch("data/reference_smn1_masked/.index.done"),
    log: "logs/smn1_mask_index.log"
    threads: 8
    resources:
        mem_mb=32000, runtime=240,
    shell:
        """
        bash scripts/align_01_mask_index.sh 2>&1 | tee {log}
        """

rule bismark_align_masked:
    input:
        r1  = "data/processed/{sample}_1_val_1.fq.gz",
        r2  = "data/processed/{sample}_2_val_2.fq.gz",
        idx = "data/reference_smn1_masked/.index.done",
    output:
        bam  = "results/alignments_smn1_masked/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        done = touch("results/alignments_smn1_masked/bs/.{sample}.align.done"),
    log: "logs/smn1_masked_align_{sample}.log"
    threads: 24
    resources:
        mem_mb=200000, runtime=2880, constraint="ehc",
    shell:
        """
        bismark --genome data/reference_smn1_masked \
            --parallel 8 -p 2 -N 1 -L 20 --score-min L,0,-0.6 \
            -1 {input.r1} -2 {input.r2} \
            -o results/alignments_smn1_masked/bs/ \
            --temp_dir results/alignments_smn1_masked/bs/tmp_{wildcards.sample} \
            2>&1 | tee {log}
        """

rule dedup_and_extract_masked:
    input:
        bam  = "results/alignments_smn1_masked/bs/{sample}_1_val_1_bismark_bt2_pe.bam",
        done = "results/alignments_smn1_masked/bs/.{sample}.align.done",
    output:
        cx   = "results/alignments_smn1_masked/cx_report/{sample}.CX_report.txt.gz",
        done = touch("results/alignments_smn1_masked/cx_report/.{sample}.cx.done"),
    log: "logs/smn1_masked_dedup_{sample}.log"
    threads: 8
    resources:
        mem_mb=64000, runtime=1440,
    shell:
        """
        mkdir -p results/alignments_smn1_masked/dedup \
                 results/alignments_smn1_masked/methylation
        deduplicate_bismark --paired --bam \
            --output_dir results/alignments_smn1_masked/dedup/ \
            {input.bam} 2>&1 | tee {log}
        DEDUP=results/alignments_smn1_masked/dedup/{wildcards.sample}_1_val_1_bismark_bt2_pe.deduplicated.bam
        bismark_methylation_extractor \
            --paired-end --CX --cytosine_report --parallel 8 \
            --genome_folder data/reference_smn1_masked \
            --output results/alignments_smn1_masked/methylation/ \
            "$DEDUP" 2>&1 | tee -a {log}
        CX_SRC=$(ls results/alignments_smn1_masked/methylation/{wildcards.sample}*CX_report.txt 2>/dev/null | head -1)
        gzip -c "$CX_SRC" > {output.cx}
        rm -f results/alignments_smn1_masked/methylation/CHG_context_{wildcards.sample}*.txt
        rm -f results/alignments_smn1_masked/methylation/CHH_context_{wildcards.sample}*.txt
        rm -f "$CX_SRC" "$DEDUP"
        """

rule split_chr5_masked:
    input:
        cx   = "results/alignments_smn1_masked/cx_report/{sample}.CX_report.txt.gz",
        done = "results/alignments_smn1_masked/cx_report/.{sample}.cx.done",
    output:
        chr5 = "results/alignments_smn1_masked/chr5_cx/{sample}_chr5.CX_report.txt",
    log: "logs/split_chr5_{sample}.log"
    shell:
        """
        mkdir -p results/alignments_smn1_masked/chr5_cx
        zcat {input.cx} | awk '$1=="chr5"' > {output.chr5} 2>&1 | tee {log}
        """

rule dmr_calling:
    input:
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
        "results/dmr/dmr_summary.tsv",
    log: "logs/dmr_calling.log"
    threads: 4
    resources:
        mem_mb=64000, runtime=360,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/02_dmr_calling.R 2>&1 | tee {log}
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/03_dmr_combine.R 2>&1 | tee -a {log}
        """

rule qc:
    input:
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        "results/dmr_qc/sample_PCA_12samples_chr1.pdf",
        "results/dmr_qc/per_sample_methylation_violin_chr1.pdf",
        "results/qc/coverage_4lines/coverage_4lines_per_condition.pdf",
        "results/qc/additional_qc/mbias_CpG_all_samples.pdf",
        "results/qc/additional_qc/bisulfite_conversion_efficiency.pdf",
    log: "logs/qc.log"
    resources:
        mem_mb=64000, runtime=240,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/01_qc.R 2>&1 | tee {log}
        """

rule dmr_annotate:
    input:
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
    output:
        expand("results/dmr_annotation/{contrast}_annotated.csv", contrast=CONTRASTS),
        expand("results/dmr_annotation/{contrast}_GO_BP_hypo_dotplot.pdf", contrast=CONTRASTS),
        "results/dmr_annotation/DMR_annotation_combined_count.pdf",
        "results/dmr_overlap/dmr_upset_plot_5contrasts.pdf",
    log: "logs/dmr_annotate.log"
    resources:
        mem_mb=32000, runtime=180,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/02_dmr_annotate.R 2>&1 | tee {log}
        """

rule enrichment:
    input:
        expand("results/dmr_annotation/{contrast}_annotated.csv", contrast=CONTRASTS),
    output:
        "results/dmr_annotation/msigdb_v2/ASO_specific_msigdb_all_combined.pdf",
        "results/h3k9me2_overlap/h3k9me2_signal_boxplot.pdf",
        "results/smn2_enhancer/smn2_h3k27ac_peak_summary.csv",
    log: "logs/enrichment.log"
    resources:
        mem_mb=32000, runtime=120,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/03_enrichment.R 2>&1 | tee {log}
        """

rule tf_motif:
    input:
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
    output:
        "results/tf_motif/motif_enrichment_volcano.pdf",
        "results/tf_motif/motif_enrichment_results.rds",
    log: "logs/tf_motif.log"
    resources:
        mem_mb=32000, runtime=120,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/04_tf_motif.R 2>&1 | tee {log}
        """

rule locus_plots:
    input:
        "results/dmr/meth_pooled_cache.rds",
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
    output:
        expand("results/dmr/plots/annotated/{gene}_annotated.pdf",
               gene=["MTA1-DT","SLC32A1","CHRNB3","GLRA4","GFRA2",
                     "SEMA3C","PHACTR3","SOX5","RNF169","SMN2"]),
        "results/lowres_profiles/lowres_allgroups_chrX_500kb.pdf",
        "results/lowres_profiles/lowres_allgroups_chr5_SMN2_10kb.pdf",
    log: "logs/locus_plots.log"
    resources:
        mem_mb=64000, runtime=180,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/05_locus_plots.R 2>&1 | tee {log}
        """

rule smn2_locus:
    input:
        expand("results/alignments_smn1_masked/chr5_cx/{sample}_chr5.CX_report.txt", sample=SAMPLES),
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        "results/smn2_locus_final/SMN_locus_masked_all_comparisons.pdf",
        "results/smn2_locus_final/SMN2_masked_vs_unmasked_comparison.pdf",
    log: "logs/smn2_locus.log"
    resources:
        mem_mb=64000, runtime=180,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/05_smn2_locus_final.R 2>&1 | tee {log}
        """

rule smn2_sensitive:
    input:
        expand("results/alignments_smn1_masked/chr5_cx/{sample}_chr5.CX_report.txt", sample=SAMPLES),
        expand("results/dmr_annotation/{contrast}_annotated.csv", contrast=CONTRASTS),
    output:
        "results/smn2_local_dmr/SMN2_sensitive_DMRs_all_contrasts.csv",
        "results/dmr_annotation/top20_hypo_ASO_by_methdiff.csv",
    log: "logs/smn2_sensitive.log"
    resources:
        mem_mb=32000, runtime=120,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/07_smn2_sensitive.R 2>&1 | tee {log}
        """

rule dss_validation:
    input:
        expand("results/alignments/bs/by_chr/.{sample}.split.done", sample=SAMPLES),
    output:
        "results/dss_replicate/DSS_GO_dotplot.pdf",
    log: "logs/dss_validation.log"
    resources:
        mem_mb=64000, runtime=240,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/08_dss_validation.R 2>&1 | tee {log}
        """

rule manhattan:
    input:
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
    output:
        expand("results/manhattan/manhattan_{contrast}.pdf", contrast=CONTRASTS),
    log: "logs/manhattan.log"
    resources:
        mem_mb=32000, runtime=60,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/09_manhattan.R 2>&1 | tee {log}
        """

rule tss_heatmap:
    input:
        "results/dmr/meth_pooled_cache.rds",
        expand("results/dmr/dmr_{contrast}.rds", contrast=CONTRASTS),
    output:
        "results/tss_metaplot/TSS_metaplot.pdf",
        "results/tss_metaplot/DMR_heatmap_top500_ASO_VPA_methdiff.pdf",
    log: "logs/tss_heatmap.log"
    resources:
        mem_mb=64000, runtime=180,
    shell:
        """
        R_LIBS_USER={R_LIBS} {RSCRIPT} scripts/06_tss_heatmap.R 2>&1 | tee {log}
        """
