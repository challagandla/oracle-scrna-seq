"""
Rule: Cell Type Annotation
"""

rule annotate:
    input:
        h5ad = "results/merged/07_markers/with_markers.h5ad"
    output:
        h5ad = "results/merged/08_annotate/annotated.h5ad",
        figs = directory("results/merged/08_annotate/figures")
    params:
        markers_json     = config["annotation"]["markers_json"] or "",
        celltypist       = "--celltypist" if config["annotation"]["celltypist"] else "",
        celltypist_model = config["annotation"]["celltypist_model"],
        manual_tsv       = config["annotation"]["manual_tsv"] or "",
    log:
        "logs/merged/08_annotate.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/08_annotate.py \
            {input.h5ad} \
            --out $(dirname {output.h5ad}) \
            {params.celltypist} \
            --celltypist-model {params.celltypist_model} \
            $([ -n "{params.markers_json}" ] && echo "--markers-json {params.markers_json}") \
            $([ -n "{params.manual_tsv}" ]   && echo "--manual-tsv {params.manual_tsv}") \
        > {log} 2>&1
        """
