"""
Rule: Quality Control — per sample
"""

rule qc:
    input:
        h5 = get_sample_path
    output:
        h5ad  = "results/{sample}/00_qc/qc_filtered.h5ad",
        stats = "results/{sample}/00_qc/qc_stats.json",
        figs  = directory("results/{sample}/00_qc/figures")
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
        "logs/{sample}/00_qc.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/00_qc.py \
            {input.h5} \
            --out $(dirname {output.h5ad}) \
            --mt-prefix {params.mt_prefix} \
            --ribo-prefix {params.ribo_prefix} \
            --hb-pattern "{params.hb_pattern}" \
            --mad-counts {params.mad_counts} \
            --mad-genes {params.mad_genes} \
            --mad-mt {params.mad_mt} \
            --mt-hard {params.mt_hard} \
            --min-cells {params.min_cells} \
            --min-counts {params.min_counts} \
            --min-genes {params.min_genes} \
        > {log} 2>&1
        """
