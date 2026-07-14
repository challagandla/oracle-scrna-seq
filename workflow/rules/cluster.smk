"""
Rule: Leiden Clustering
"""

rule cluster:
    input:
        h5ad = result_path("merged/05_embedding/embedded.h5ad"),
        code = [
            "scripts/06_cluster.py",
            "scripts/utils/io_utils.py",
            "scripts/utils/plot_utils.py",
        ]
    output:
        h5ad    = result_path("merged/06_cluster/clustered.h5ad"),
        summary = result_path("merged/06_cluster/clustering_summary.csv"),
        figs    = directory(result_path("merged/06_cluster/figures"))
    threads: resource_value("cluster", "threads")
    resources:
        mem_mb=resource_value("cluster", "mem_mb")
    params:
        resolutions = " ".join(str(r) for r in config["clustering"]["resolutions"]),
        default_res = config["clustering"]["default_resolution"],
        use_rep     = use_rep(),
        seed        = config.get("random_seed", 42),
    log:
        log_path("merged/06_cluster.log")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/06_cluster.py \
            {input.h5ad:q} \
            --out "$(dirname {output.h5ad:q})" \
            --resolutions {params.resolutions} \
            --default-resolution {params.default_res} \
            --use-rep {params.use_rep:q} \
            --seed {params.seed} \
        > {log:q} 2>&1
        """
