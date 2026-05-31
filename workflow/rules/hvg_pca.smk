"""
Rule: HVG Selection & PCA
"""

rule hvg_pca:
    input:
        h5ad = "results/merged/02_normalize/normalized.h5ad"
    output:
        h5ad = "results/merged/03_hvg_pca/hvg_pca.h5ad",
        figs = directory("results/merged/03_hvg_pca/figures")
    params:
        n_top_genes = config["hvg_pca"]["n_top_genes"],
        flavor      = config["hvg_pca"]["flavor"],
        n_pcs       = config["hvg_pca"]["n_pcs"],
        batch_key   = config["batch"]["batch_key"] or "",
    log:
        "logs/merged/03_hvg_pca.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/03_hvg_pca.py \
            {input.h5ad} \
            --out $(dirname {output.h5ad}) \
            --n-hvg {params.n_top_genes} \
            --hvg-flavor {params.flavor} \
            --n-pcs {params.n_pcs} \
            $([ -n "{params.batch_key}" ] && echo "--batch-key {params.batch_key}") \
        > {log} 2>&1
        """
