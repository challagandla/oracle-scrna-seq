"""
Rules: optional raw FASTQ QC, FastQ Screen contamination detection, and MultiQC.

The scRNA workflow can start from processed count matrices only. These rules run
only for samples that provide raw_fastq_r1/raw_fastq_r2 in config/samples.tsv.
"""

rule fastqc_raw_fastq:
    input:
        R1=lambda wc: get_raw_fastq(wc, "raw_fastq_r1"),
        R2=lambda wc: get_raw_fastq(wc, "raw_fastq_r2")
    output:
        R1html="results/{sample}/raw_fastq_qc/fastqc/{sample}_R1_fastqc.html",
        R1zip="results/{sample}/raw_fastq_qc/fastqc/{sample}_R1_fastqc.zip",
        R2html="results/{sample}/raw_fastq_qc/fastqc/{sample}_R2_fastqc.html",
        R2zip="results/{sample}/raw_fastq_qc/fastqc/{sample}_R2_fastqc.zip"
    threads: 2
    log:
        "logs/{sample}/raw_fastq_qc/fastqc.log"
    conda:
        "../../envs/raw_fastq_qc.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}/raw_fastq_qc/fastqc logs/{wildcards.sample}/raw_fastq_qc
        tmpdir=$(mktemp -d results/{wildcards.sample}/raw_fastq_qc/fastqc/tmp.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        ln -sf "$(readlink -f {input.R1})" "$tmpdir/{wildcards.sample}_R1.fastq.gz"
        ln -sf "$(readlink -f {input.R2})" "$tmpdir/{wildcards.sample}_R2.fastq.gz"
        fastqc \
            -o results/{wildcards.sample}/raw_fastq_qc/fastqc \
            -t {threads} \
            "$tmpdir/{wildcards.sample}_R1.fastq.gz" \
            "$tmpdir/{wildcards.sample}_R2.fastq.gz" \
            > {log} 2>&1
        """


rule fastq_screen_raw_fastq:
    input:
        R1=lambda wc: get_raw_fastq(wc, "raw_fastq_r1"),
        R2=lambda wc: get_raw_fastq(wc, "raw_fastq_r2"),
        conf=FASTQ_SCREEN_CONF
    output:
        R1txt="results/{sample}/raw_fastq_qc/fastq_screen/{sample}_R1_screen.txt",
        R1html="results/{sample}/raw_fastq_qc/fastq_screen/{sample}_R1_screen.html",
        R2txt="results/{sample}/raw_fastq_qc/fastq_screen/{sample}_R2_screen.txt",
        R2html="results/{sample}/raw_fastq_qc/fastq_screen/{sample}_R2_screen.html"
    params:
        subset=FASTQ_SCREEN_SUBSET
    threads: 4
    log:
        "logs/{sample}/raw_fastq_qc/fastq_screen.log"
    conda:
        "../../envs/raw_fastq_qc.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}/raw_fastq_qc/fastq_screen logs/{wildcards.sample}/raw_fastq_qc
        tmpdir=$(mktemp -d results/{wildcards.sample}/raw_fastq_qc/fastq_screen/tmp.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        ln -sf "$(readlink -f {input.R1})" "$tmpdir/{wildcards.sample}_R1.fastq.gz"
        ln -sf "$(readlink -f {input.R2})" "$tmpdir/{wildcards.sample}_R2.fastq.gz"
        fastq_screen \
            --conf {input.conf} \
            --aligner bowtie2 \
            --threads {threads} \
            --subset {params.subset} \
            --outdir results/{wildcards.sample}/raw_fastq_qc/fastq_screen \
            "$tmpdir/{wildcards.sample}_R1.fastq.gz" \
            "$tmpdir/{wildcards.sample}_R2.fastq.gz" \
            > {log} 2>&1
        """


rule multiqc:
    input:
        raw_qc=RAW_FASTQC_ZIP + FASTQ_SCREEN_TEXT + FASTQ_SCREEN_HTML,
        html_report="results/report.html"
    output:
        html=MULTIQC_REPORT
    params:
        run_multiqc="1" if RAW_FASTQ_SAMPLES else "0"
    log:
        "logs/multiqc.log"
    conda:
        "../../envs/raw_fastq_qc.yaml"
    shell:
        """
        mkdir -p results/multiqc logs
        if [ "{params.run_multiqc}" = "1" ]; then
            multiqc results logs --outdir results/multiqc --filename multiqc_report.html --force > {log} 2>&1 || \
            printf '%s\n' '<html><body><h1>MultiQC</h1><p>MultiQC did not find supported raw FASTQ QC modules. Check logs/multiqc.log.</p></body></html>' > {output.html}
        else
            printf '%s\n' '<html><body><h1>MultiQC</h1><p>No raw FASTQ files were configured in config/samples.tsv, so FastQC/FastQ Screen were skipped.</p></body></html>' > {output.html}
            printf '%s\n' 'No raw FASTQ files were configured; wrote placeholder MultiQC report.' > {log}
        fi
        """
