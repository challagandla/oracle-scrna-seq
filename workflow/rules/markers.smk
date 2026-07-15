"""
Rule: Marker Genes & Differential Expression
"""

rule markers:
    input:
        h5ad = result_path("merged/06_cluster/clustered.h5ad"),
        code = [
            "scripts/07_markers_de.py",
            "scripts/utils/io_utils.py",
            "scripts/utils/plot_utils.py",
        ]
    output:
        h5ad = result_path("merged/07_markers/with_markers.h5ad"),
        csv  = result_path("merged/07_markers/cluster_markers.csv"),
        status = result_path("merged/07_markers/pseudobulk_status.csv"),
        de = directory(result_path("merged/07_markers/DE")),
        figs = directory(result_path("merged/07_markers/figures"))
    threads: resource_value("markers", "threads")
    resources:
        mem_mb=resource_value("markers", "mem_mb")
    params:
        method        = config["markers"]["method"],
        n_genes       = config["markers"]["n_genes"],
        logfc_min     = config["markers"]["logfc_min"],
        pct_min       = config["markers"]["pct_min"],
        condition_options = marker_condition_options(),
        covariate_options = marker_covariate_options(),
        min_cells     = config["markers"].get("min_cells_per_sample", 10),
        min_replicates = config["markers"].get("min_replicates_per_group", 2),
        min_total_count = config["markers"].get("min_total_count", 10),
    log:
        log_path("merged/07_markers.log")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/07_markers_de.py \
            {input.h5ad:q} \
            --out "$(dirname {output.h5ad:q})" \
            --method {params.method:q} \
            --n-genes {params.n_genes} \
            --logfc-min {params.logfc_min} \
            --pct-min {params.pct_min} \
            {params.condition_options:q} \
            {params.covariate_options:q} \
            --min-cells-per-sample {params.min_cells} \
            --min-replicates-per-group {params.min_replicates} \
            --min-total-count {params.min_total_count} \
        > {log:q} 2>&1
        """
