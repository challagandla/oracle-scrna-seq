"""Optional raw-read diagnostics; quantification remains an upstream boundary."""


rule fastqc_raw_fastq:
    input:
        R1=lambda wc: get_raw_fastq(wc, "raw_fastq_r1"),
        R2=lambda wc: get_raw_fastq(wc, "raw_fastq_r2")
    output:
        R1html=result_path("{sample}/raw_fastq_qc/fastqc/{sample}_R1_fastqc.html"),
        R1zip=result_path("{sample}/raw_fastq_qc/fastqc/{sample}_R1_fastqc.zip"),
        R2html=result_path("{sample}/raw_fastq_qc/fastqc/{sample}_R2_fastqc.html"),
        R2zip=result_path("{sample}/raw_fastq_qc/fastqc/{sample}_R2_fastqc.zip")
    threads: resource_value("raw_fastqc", "threads")
    resources:
        mem_mb=resource_value("raw_fastqc", "mem_mb")
    log:
        log_path("{sample}/raw_fastq_qc/fastqc.log")
    conda:
        "../../envs/raw_fastq_qc.yaml"
    shell:
        """
        outdir=$(dirname {output.R1html:q})
        mkdir -p "$outdir" "$(dirname {log:q})"
        tmpdir=$(mktemp -d "$outdir/tmp.XXXXXX")
        trap 'rm -rf "$tmpdir"' EXIT
        r1=$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' {input.R1:q})
        r2=$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' {input.R2:q})
        ln -s "$r1" "$tmpdir/{wildcards.sample}_R1.fastq.gz"
        ln -s "$r2" "$tmpdir/{wildcards.sample}_R2.fastq.gz"
        fastqc -o "$outdir" -t {threads} \
            "$tmpdir/{wildcards.sample}_R1.fastq.gz" \
            "$tmpdir/{wildcards.sample}_R2.fastq.gz" \
            > {log:q} 2>&1
        """


rule fastq_screen_raw_fastq:
    input:
        fq=get_screen_fastq,
        conf=FASTQ_SCREEN_CONF
    output:
        txt=result_path("{sample}/raw_fastq_qc/fastq_screen/{sample}_{read}_screen.txt"),
        html=result_path("{sample}/raw_fastq_qc/fastq_screen/{sample}_{read}_screen.html")
    params:
        subset=FASTQ_SCREEN_SUBSET
    threads: resource_value("fastq_screen", "threads")
    resources:
        mem_mb=resource_value("fastq_screen", "mem_mb")
    wildcard_constraints:
        read="R[12]"
    log:
        log_path("{sample}/raw_fastq_qc/fastq_screen_{read}.log")
    conda:
        "../../envs/raw_fastq_qc.yaml"
    shell:
        """
        outdir=$(dirname {output.txt:q})
        mkdir -p "$outdir" "$(dirname {log:q})"
        tmpdir=$(mktemp -d "$outdir/tmp.XXXXXX")
        trap 'rm -rf "$tmpdir"' EXIT
        fq=$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' {input.fq:q})
        ln -s "$fq" "$tmpdir/{wildcards.sample}_{wildcards.read}.fastq.gz"
        fastq_screen \
            --conf {input.conf:q} \
            --aligner bowtie2 \
            --threads {threads} \
            --subset {params.subset} \
            --outdir "$outdir" \
            "$tmpdir/{wildcards.sample}_{wildcards.read}.fastq.gz" \
            > {log:q} 2>&1
        """


if RAW_FASTQ_SAMPLES:
    rule multiqc:
        input:
            raw_qc=RAW_FASTQC_ZIP + FASTQ_SCREEN_TEXT + FASTQ_SCREEN_HTML,
            html_report=result_path("report.html")
        output:
            html=MULTIQC_REPORT
        threads: resource_value("multiqc", "threads")
        resources:
            mem_mb=resource_value("multiqc", "mem_mb")
        log:
            log_path("multiqc.log")
        conda:
            "../../envs/raw_fastq_qc.yaml"
        shell:
            """
            outdir=$(dirname {output.html:q})
            mkdir -p "$outdir" "$(dirname {log:q})"
            multiqc {RESULTS_DIR:q} {LOGS_DIR:q} \
                --outdir "$outdir" \
                --filename multiqc_report.html \
                --force > {log:q} 2>&1
            """
else:
    rule multiqc_placeholder:
        input:
            html_report=result_path("report.html")
        output:
            html=MULTIQC_REPORT
        log:
            log_path("multiqc.log")
        run:
            Path(output.html).parent.mkdir(parents=True, exist_ok=True)
            Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
            Path(output.html).write_text(
                "<html><body><h1>Raw FASTQ QC not requested</h1>"
                "<p>No complete raw FASTQ pairs were configured. This does not "
                "assess alignment, quantification, cell calling, or ambient RNA.</p>"
                "</body></html>",
                encoding="utf-8",
            )
            Path(log[0]).write_text("No raw FASTQ pairs configured.\n", encoding="utf-8")
