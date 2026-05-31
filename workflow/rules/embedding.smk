"""
Rule: Neighbor Graph & UMAP
"""

rule embedding:
    input:
        h5ad = "results/merged/04_batch/batch_corrected.h5ad"
    output:
        h5ad = "results/merged/05_embedding/embedded.h5ad",
        figs = directory("results/merged/05_embedding/figures")
    params:
        n_neighbors = config["embedding"]["n_neighbors"],
        n_pcs       = config["hvg_pca"]["n_pcs_use"],
        min_dist    = config["embedding"]["min_dist"],
        spread      = config["embedding"]["spread"],
        use_rep     = use_rep(),
        batch_key   = config["batch"]["batch_key"] or "",
    log:
        "logs/merged/05_embedding.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/05_neighbors_umap.py \
            {input.h5ad} \
            --out $(dirname {output.h5ad}) \
            --n-neighbors {params.n_neighbors} \
            --n-pcs {params.n_pcs} \
            --min-dist {params.min_dist} \
            --spread {params.spread} \
            --use-rep {params.use_rep} \
            $([ -n "{params.batch_key}" ] && echo "--batch-key {params.batch_key}") \
        > {log} 2>&1
        """
