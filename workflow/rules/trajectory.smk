"""
Rule: Trajectory analysis (PAGA + scVelo) — optional
"""

rule trajectory:
    input:
        h5ad = result_path("merged/08_annotate/annotated.h5ad"),
        auxiliary=lambda wc: [
            path for path in (
                config["trajectory"].get("loom_file", ""),
                config["trajectory"].get("barcode_map_tsv", ""),
            ) if _filled(path)
        ],
        code = ["scripts/09_trajectory.py", "scripts/utils/io_utils.py"]
    output:
        h5ad = TRAJECTORY_H5AD,
        figs = directory(result_path("merged/09_trajectory/figures"))
    params:
        paga_option = cli_switch("--paga", config["trajectory"]["run_paga"]),
        loom_option = cli_option("--loom", config["trajectory"].get("loom_file", "")),
        barcode_map_option = cli_option(
            "--barcode-map", config["trajectory"].get("barcode_map_tsv", "")
        ),
        min_shared = config["trajectory"].get("min_shared_cells", 100),
        seed       = config.get("random_seed", 42),
    log:
        log_path("merged/09_trajectory.log")
    threads: resource_value("trajectory", "threads")
    resources:
        mem_mb=resource_value("trajectory", "mem_mb")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/09_trajectory.py \
            {input.h5ad:q} \
            --out "$(dirname {output.h5ad:q})" \
            {params.paga_option:q} \
            {params.loom_option:q} \
            {params.barcode_map_option:q} \
            --min-shared-cells {params.min_shared} \
            --threads {threads} \
            --seed {params.seed} \
        > {log:q} 2>&1
        """
