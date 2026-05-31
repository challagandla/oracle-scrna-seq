"""
Rule: Leiden Clustering
"""

rule cluster:
    input:
        h5ad = "results/merged/05_embedding/embedded.h5ad"
    output:
        h5ad    = "results/merged/06_cluster/clustered.h5ad",
        summary = "results/merged/06_cluster/clustering_summary.csv",
        figs    = directory("results/merged/06_cluster/figures")
    params:
        resolutions = " ".join(str(r) for r in config["clustering"]["resolutions"]),
        default_res = config["clustering"]["default_resolution"],
        use_rep     = use_rep(),
    log:
        "logs/merged/06_cluster.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/06_cluster.py \
            {input.h5ad} \
            --out $(dirname {output.h5ad}) \
            --resolutions {params.resolutions} \
            --default-resolution {params.default_res} \
            --use-rep {params.use_rep} \
        > {log} 2>&1
        """
