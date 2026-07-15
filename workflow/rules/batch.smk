"""
Rule: Batch Correction (skipped if batch_key is empty)
"""

rule batch_correct:
    input:
        h5ad = result_path("merged/03_hvg_pca/hvg_pca.h5ad"),
        code = [
            "scripts/04_batch_correct.py",
            "scripts/utils/io_utils.py",
            "scripts/utils/plot_utils.py",
        ]
    output:
        h5ad = result_path("merged/04_batch/batch_corrected.h5ad"),
        figs = directory(result_path("merged/04_batch/figures"))
    threads: resource_value("batch", "threads")
    resources:
        mem_mb=resource_value("batch", "mem_mb")
    params:
        batch_option = cli_option("--batch-key", config["batch"]["batch_key"]),
        method    = config["batch"]["method"],
        n_latent  = config["batch"]["n_latent"],
        seed      = config.get("random_seed", 42),
    log:
        log_path("merged/04_batch.log")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/04_batch_correct.py \
            {input.h5ad:q} \
            --out "$(dirname {output.h5ad:q})" \
            --method {params.method:q} \
            --n-latent {params.n_latent} \
            {params.batch_option:q} \
            --seed {params.seed} \
        > {log:q} 2>&1
        """
