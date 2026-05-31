"""
Rule: HTML summary report
"""

rule report:
    input:
        annotated = "results/merged/08_annotate/annotated.h5ad",
        qc_figs   = expand("results/{sample}/00_qc/figures", sample=SAMPLES),
        clust_csv = "results/merged/06_cluster/clustering_summary.csv",
        markers   = "results/merged/07_markers/cluster_markers.csv",
    output:
        html = "results/report.html"
    log:
        "logs/report.log"
    conda:
        "../../envs/scrna.yaml"
    shell:
        """
        python scripts/make_report.py \
            --annotated {input.annotated} \
            --markers   {input.markers} \
            --clust-csv {input.clust_csv} \
            --samples-tsv config/samples.tsv \
            --out {output.html} \
        > {log} 2>&1
        """
