"""
Rule: Doublet Detection — per sample
"""

rule doublets:
    input:
        h5ad = "results/{sample}/00_qc/qc_filtered.h5ad"
    output:
        h5ad = "results/{sample}/01_doublets/no_doublets.h5ad",
        figs = directory("results/{sample}/01_doublets/figures")
    params:
        method        = config["doublets"]["method"],
        expected_rate = config["doublets"]["expected_rate"],
    log:
        "logs/{sample}/01_doublets.log"
    conda:
        lambda wildcards: (
            "../../envs/r_env.yaml"
            if config["doublets"]["method"] == "scdblfinder"
            else "../../envs/scrna.yaml"
        )
    shell:
        """
        python scripts/01_doublets.py \
            {input.h5ad} \
            --out $(dirname {output.h5ad}) \
            --method {params.method} \
            --expected-rate {params.expected_rate} \
        > {log} 2>&1
        """
