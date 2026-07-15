"""
Rule: HTML summary report
"""

rule snapshot_effective_config:
    output:
        effective = result_path("run_config.yaml")
    params:
        config_digest = CONFIG_DIGEST
    run:
        Path(output.effective).parent.mkdir(parents=True, exist_ok=True)
        Path(output.effective).write_text(
            yaml.safe_dump(dict(config), sort_keys=False),
            encoding="utf-8",
        )

rule report:
    input:
        annotated = result_path("merged/08_annotate/annotated.h5ad"),
        qc_figs   = expand(result_path("{sample}/00_qc/figures"), sample=SAMPLES),
        clust_csv = result_path("merged/06_cluster/clustering_summary.csv"),
        markers   = result_path("merged/07_markers/cluster_markers.csv"),
        de_status = result_path("merged/07_markers/pseudobulk_status.csv"),
        annotation = result_path("merged/08_annotate/annotation_summary.csv"),
        marker_coverage = result_path("merged/08_annotate/marker_coverage.csv"),
        samples_tsv = config["samples_tsv"],
        config_yaml = result_path("run_config.yaml"),
        environment_specs = [
            "envs/scrna.yaml",
            "envs/r_env.yaml",
            "envs/raw_fastq_qc.yaml",
        ],
        downstream = TRAJECTORY_TARGETS,
        code = ["scripts/make_report.py"]
    output:
        html = result_path("report.html"),
        manifest = result_path("run_manifest.json")
    threads: resource_value("report", "threads")
    resources:
        mem_mb=resource_value("report", "mem_mb")
    params:
        results_dir=RESULTS_DIR,
        hash_inputs="--hash-inputs" if config.get("provenance", {}).get("hash_inputs", False) else "",
        git_commit=GIT_COMMIT
    log:
        log_path("report.log")
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        python scripts/make_report.py \
            --annotated {input.annotated:q} \
            --markers {input.markers:q} \
            --clust-csv {input.clust_csv:q} \
            --annotation-summary {input.annotation:q} \
            --marker-coverage {input.marker_coverage:q} \
            --de-status {input.de_status:q} \
            --samples-tsv {input.samples_tsv:q} \
            --config-yaml {input.config_yaml:q} \
            --environment-specs {input.environment_specs:q} \
            --results-dir {params.results_dir:q} \
            --out {output.html:q} \
            --manifest {output.manifest:q} \
            {params.hash_inputs} \
        > {log:q} 2>&1
        """
