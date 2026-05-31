"""
Rule: Trajectory analysis (PAGA + scVelo) — optional
"""

rule trajectory:
    input:
        h5ad = "results/merged/08_annotate/annotated.h5ad"
    output:
        h5ad = "results/merged/09_trajectory/trajectory.h5ad",
        figs = directory("results/merged/09_trajectory/figures")
    params:
        paga      = "--paga" if config["trajectory"]["run_paga"] else "",
        loom_file = config["trajectory"]["loom_file"] or "",
    log:
        "logs/merged/09_trajectory.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/09_trajectory.py \
            {input.h5ad} \
            --out $(dirname {output.h5ad}) \
            {params.paga} \
            $([ -n "{params.loom_file}" ] && echo "--loom {params.loom_file}") \
        > {log} 2>&1
        """
