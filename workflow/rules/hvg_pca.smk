"""
Rule: HVG Selection & PCA
"""

rule hvg_pca:
    input:
        h5ad = result_path("merged/02_normalize/normalized.h5ad"),
        code = [
            "scripts/03_hvg_pca.py",
            "scripts/utils/io_utils.py",
            "scripts/utils/plot_utils.py",
        ]
    output:
        h5ad = result_path("merged/03_hvg_pca/hvg_pca.h5ad"),
        figs = directory(result_path("merged/03_hvg_pca/figures"))
    threads: resource_value("hvg_pca", "threads")
    resources:
        mem_mb=resource_value("hvg_pca", "mem_mb")
    params:
        n_top_genes = config["hvg_pca"]["n_top_genes"],
        flavor      = config["hvg_pca"]["flavor"],
        n_pcs       = config["hvg_pca"]["n_pcs"],
        n_pcs_use   = config["hvg_pca"]["n_pcs_use"],
        batch_option = cli_option("--batch-key", config["batch"]["batch_key"]),
        seed        = config.get("random_seed", 42),
    log:
        log_path("merged/03_hvg_pca.log")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/03_hvg_pca.py \
            {input.h5ad:q} \
            --out "$(dirname {output.h5ad:q})" \
            --n-hvg {params.n_top_genes} \
            --hvg-flavor {params.flavor} \
            --n-pcs {params.n_pcs} \
            --n-pcs-use {params.n_pcs_use} \
            {params.batch_option:q} \
            --seed {params.seed} \
        > {log:q} 2>&1
        """
