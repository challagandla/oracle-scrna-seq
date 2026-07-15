"""
Rule: Doublet Detection — per sample
"""

rule doublets:
    input:
        h5ad = result_path("{sample}/00_qc/qc_filtered.h5ad"),
        code = [
            "scripts/01_doublets.py",
            "scripts/utils/io_utils.py",
            "scripts/utils/validation.py",
        ]
    output:
        h5ad = result_path("{sample}/01_doublets/no_doublets.h5ad"),
        figs = directory(result_path("{sample}/01_doublets/figures"))
    threads: resource_value("doublets", "threads")
    resources:
        mem_mb=resource_value("doublets", "mem_mb")
    params:
        method        = config["doublets"]["method"],
        expected_rate = config["doublets"]["expected_rate"],
        seed          = config.get("random_seed", 42),
    log:
        log_path("{sample}/01_doublets.log")
    conda:
        lambda wildcards: (
            "../../envs/r_env.yaml"
            if config["doublets"]["method"] == "scdblfinder"
            else "../../envs/scrna.yaml"
        )
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/01_doublets.py \
            {input.h5ad:q} \
            --out "$(dirname {output.h5ad:q})" \
            --method {params.method:q} \
            --expected-rate {params.expected_rate} \
            --seed {params.seed} \
        > {log:q} 2>&1
        """
