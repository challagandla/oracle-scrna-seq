"""
Rule: Merge samples → Normalize
"""

rule merge_samples:
    input:
        expand("results/{sample}/01_doublets/no_doublets.h5ad", sample=SAMPLES)
    output:
        h5ad = "results/merged/00_merged/merged.h5ad"
    params:
        samples_tsv = config["samples_tsv"]
    log:
        "logs/merged/00_merge.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/merge_samples.py \
            --samples-tsv {params.samples_tsv} \
            --input-dir results \
            --out $(dirname {output.h5ad}) \
        > {log} 2>&1
        """

rule normalize:
    input:
        h5ad = "results/merged/00_merged/merged.h5ad"
    output:
        h5ad = "results/merged/02_normalize/normalized.h5ad",
        figs = directory("results/merged/02_normalize/figures")
    params:
        method     = config["normalization"]["method"],
        target_sum = config["normalization"]["target_sum"],
    log:
        "logs/merged/02_normalize.log"
    conda:
        lambda wildcards: (
            "../../envs/r_env.yaml"
            if config["normalization"]["method"] == "scran"
            else "../../envs/scrna.yaml"
        )
    shell:
        """
        python scripts/02_normalize.py \
            {input.h5ad} \
            --out $(dirname {output.h5ad}) \
            --method {params.method} \
            --target-sum {params.target_sum} \
        > {log} 2>&1
        """
