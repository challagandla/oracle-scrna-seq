"""
Rule: Quality Control — per sample
"""

rule qc:
    input:
        h5 = get_sample_path,
        code = [
            "scripts/00_qc.py",
            "scripts/utils/io_utils.py",
            "scripts/utils/plot_utils.py",
            "scripts/utils/qc_utils.py",
            "scripts/utils/validation.py",
        ]
    output:
        h5ad  = result_path("{sample}/00_qc/qc_filtered.h5ad"),
        stats = result_path("{sample}/00_qc/qc_stats.json"),
        figs  = directory(result_path("{sample}/00_qc/figures"))
    threads: resource_value("qc", "threads")
    resources:
        mem_mb=resource_value("qc", "mem_mb")
    params:
        mt_prefix   = lambda _: config.get("mt_prefix") or
                      ("mt-" if config["species"] == "mouse" else "MT-"),
        ribo_prefix = lambda _: config.get("ribo_prefix") or
                      ("Rpl,Rps" if config["species"] == "mouse" else "RPL,RPS"),
        hb_pattern  = lambda _: config.get("hb_pattern") or
                      ("^Hb[^(p)]" if config["species"] == "mouse" else "^HB[^(P)]"),
        mad_counts  = config["qc"]["mad_counts"],
        mad_genes   = config["qc"]["mad_genes"],
        mad_mt      = config["qc"]["mad_mt"],
        mt_hard     = config["qc"]["mt_hard"],
        min_cells   = config["qc"]["min_cells"],
        min_counts  = config["qc"]["min_counts"],
        min_genes   = config["qc"]["min_genes"],
    log:
        log_path("{sample}/00_qc.log")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/00_qc.py \
            {input.h5:q} \
            --out "$(dirname {output.h5ad:q})" \
            --mt-prefix {params.mt_prefix:q} \
            --ribo-prefix {params.ribo_prefix:q} \
            --hb-pattern {params.hb_pattern:q} \
            --mad-counts {params.mad_counts} \
            --mad-genes {params.mad_genes} \
            --mad-mt {params.mad_mt} \
            --mt-hard {params.mt_hard} \
            --min-cells 0 \
            --min-counts {params.min_counts} \
            --min-genes {params.min_genes} \
        > {log:q} 2>&1
        """
