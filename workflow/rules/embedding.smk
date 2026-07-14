"""
Rule: Neighbor Graph & UMAP
"""

rule embedding:
    input:
        h5ad = result_path("merged/04_batch/batch_corrected.h5ad"),
        code = [
            "scripts/05_neighbors_umap.py",
            "scripts/utils/io_utils.py",
            "scripts/utils/plot_utils.py",
        ]
    output:
        h5ad = result_path("merged/05_embedding/embedded.h5ad"),
        figs = directory(result_path("merged/05_embedding/figures"))
    threads: resource_value("embedding", "threads")
    resources:
        mem_mb=resource_value("embedding", "mem_mb")
    params:
        n_neighbors = config["embedding"]["n_neighbors"],
        n_pcs       = config["hvg_pca"]["n_pcs_use"],
        min_dist    = config["embedding"]["min_dist"],
        spread      = config["embedding"]["spread"],
        metric      = config["embedding"]["metric"],
        use_rep     = use_rep(),
        batch_option = cli_option("--batch-key", config["batch"]["batch_key"]),
        seed        = config.get("random_seed", 42),
    log:
        log_path("merged/05_embedding.log")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/05_neighbors_umap.py \
            {input.h5ad:q} \
            --out "$(dirname {output.h5ad:q})" \
            --n-neighbors {params.n_neighbors} \
            --n-pcs {params.n_pcs} \
            --min-dist {params.min_dist} \
            --spread {params.spread} \
            --metric {params.metric:q} \
            --use-rep {params.use_rep:q} \
            {params.batch_option:q} \
            --seed {params.seed} \
        > {log:q} 2>&1
        """
