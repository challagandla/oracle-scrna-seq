"""
Rule: Cell Type Annotation
"""

rule annotate:
    input:
        h5ad = result_path("merged/07_markers/with_markers.h5ad"),
        auxiliary=lambda wc: [
            path for path in (
                config["annotation"].get("markers_json", ""),
                config["annotation"].get("manual_tsv", ""),
            ) if _filled(path)
        ],
        code = [
            "scripts/08_annotate.py",
            "scripts/utils/io_utils.py",
            "scripts/utils/plot_utils.py",
        ]
    output:
        h5ad = result_path("merged/08_annotate/annotated.h5ad"),
        summary = result_path("merged/08_annotate/annotation_summary.csv"),
        coverage = result_path("merged/08_annotate/marker_coverage.csv"),
        composition = result_path("merged/08_annotate/composition_by_sample.csv"),
        figs = directory(result_path("merged/08_annotate/figures"))
    threads: resource_value("annotate", "threads")
    resources:
        mem_mb=resource_value("annotate", "mem_mb")
    params:
        marker_set_option = cli_option(
            "--marker-set", config["annotation"].get("marker_set", "")
        ),
        markers_json_option = cli_option(
            "--markers-json", config["annotation"].get("markers_json", "")
        ),
        celltypist_option = cli_switch(
            "--celltypist", config["annotation"].get("celltypist", False)
        ),
        celltypist_model = config["annotation"]["celltypist_model"],
        manual_tsv_option = cli_option(
            "--manual-tsv", config["annotation"].get("manual_tsv", "")
        ),
        min_coverage     = config["annotation"].get("min_marker_coverage", 0.5),
        min_margin       = config["annotation"].get("min_score_margin", 0.1),
        min_celltypist   = config["annotation"].get("celltypist_min_confidence", 0.5),
        seed              = config.get("random_seed", 42),
    log:
        log_path("merged/08_annotate.log")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/08_annotate.py \
            {input.h5ad:q} \
            --out "$(dirname {output.h5ad:q})" \
            {params.celltypist_option:q} \
            {params.marker_set_option:q} \
            {params.markers_json_option:q} \
            {params.manual_tsv_option:q} \
            --celltypist-model {params.celltypist_model:q} \
            --min-marker-coverage {params.min_coverage} \
            --min-score-margin {params.min_margin} \
            --celltypist-min-confidence {params.min_celltypist} \
            --seed {params.seed} \
        > {log:q} 2>&1
        """
