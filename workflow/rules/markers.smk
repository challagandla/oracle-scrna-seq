"""
Rule: Marker Genes & Differential Expression
"""

rule markers:
    input:
        h5ad = "results/merged/06_cluster/clustered.h5ad"
    output:
        h5ad = "results/merged/07_markers/with_markers.h5ad",
        csv  = "results/merged/07_markers/cluster_markers.csv",
        figs = directory("results/merged/07_markers/figures")
    params:
        method        = config["markers"]["method"],
        n_genes       = config["markers"]["n_genes"],
        logfc_min     = config["markers"]["logfc_min"],
        pct_min       = config["markers"]["pct_min"],
        condition_key = config["markers"]["condition_key"] or "",
        group1        = config["markers"]["group1"] or "",
        group2        = config["markers"]["group2"] or "",
    log:
        "logs/merged/07_markers.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/07_markers_de.py \
            {input.h5ad} \
            --out $(dirname {output.h5ad}) \
            --method {params.method} \
            --n-genes {params.n_genes} \
            --logfc-min {params.logfc_min} \
            --pct-min {params.pct_min} \
            $([ -n "{params.condition_key}" ] && echo \
              "--condition-key {params.condition_key} \
               --group1 {params.group1} \
               --group2 {params.group2}") \
        > {log} 2>&1
        """
