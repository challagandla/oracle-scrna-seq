"""
Rule: Merge samples → Normalize
"""

rule merge_samples:
    input:
        h5ads=expand(result_path("{sample}/01_doublets/no_doublets.h5ad"), sample=SAMPLES),
        samples_tsv=config["samples_tsv"],
        code=["scripts/merge_samples.py", "scripts/utils/io_utils.py"]
    output:
        h5ad = result_path("merged/00_merged/merged.h5ad")
    threads: resource_value("merge", "threads")
    resources:
        mem_mb=resource_value("merge", "mem_mb")
    params:
        input_dir=RESULTS_DIR,
        min_cells=config["qc"]["min_cells"]
    log:
        log_path("merged/00_merge.log")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/merge_samples.py \
            --samples-tsv {input.samples_tsv:q} \
            --input-dir {params.input_dir:q} \
            --min-cells {params.min_cells} \
            --out "$(dirname {output.h5ad:q})" \
        > {log:q} 2>&1
        """

rule normalize:
    input:
        h5ad = result_path("merged/00_merged/merged.h5ad"),
        code = [
            "scripts/02_normalize.py",
            "scripts/utils/io_utils.py",
            "scripts/utils/validation.py",
        ]
    output:
        h5ad = result_path("merged/02_normalize/normalized.h5ad"),
        figs = directory(result_path("merged/02_normalize/figures"))
    threads: resource_value("normalize", "threads")
    resources:
        mem_mb=resource_value("normalize", "mem_mb")
    params:
        method     = config["normalization"]["method"],
        target_sum = config["normalization"]["target_sum"],
        seed       = config.get("random_seed", 42),
    log:
        log_path("merged/02_normalize.log")
    conda:
        lambda wildcards: (
            "../../envs/r_env.yaml"
            if config["normalization"]["method"] == "scran"
            else "../../envs/scrna.yaml"
        )
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/02_normalize.py \
            {input.h5ad:q} \
            --out "$(dirname {output.h5ad:q})" \
            --method {params.method:q} \
            --target-sum {params.target_sum} \
            --seed {params.seed} \
        > {log:q} 2>&1
        """
