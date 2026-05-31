"""
Rule: Batch Correction (skipped if batch_key is empty)
"""

rule batch_correct:
    input:
        h5ad = "results/merged/03_hvg_pca/hvg_pca.h5ad"
    output:
        h5ad = "results/merged/04_batch/batch_corrected.h5ad",
        figs = directory("results/merged/04_batch/figures")
    params:
        batch_key = config["batch"]["batch_key"] or "",
        method    = config["batch"]["method"],
        n_latent  = config["batch"]["n_latent"],
    log:
        "logs/merged/04_batch.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/04_batch_correct.py \
            {input.h5ad} \
            --out $(dirname {output.h5ad}) \
            --method {params.method} \
            --n-latent {params.n_latent} \
            $([ -n "{params.batch_key}" ] && echo "--batch-key {params.batch_key}") \
        > {log} 2>&1
        """
