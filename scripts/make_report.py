#!/usr/bin/env python3
"""
Generate a self-contained HTML summary report for the pipeline run.
"""

import argparse, base64, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import pandas as pd

TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>scRNA-seq Pipeline Report</title>
  <style>
    body {{ font-family: 'Segoe UI', Arial, sans-serif; margin:0; background:#f5f5f5; color:#222; }}
    header {{ background:#2c3e50; color:#fff; padding:20px 40px; }}
    header h1 {{ margin:0; font-size:1.6em; }}
    header p  {{ margin:4px 0 0; opacity:0.7; font-size:0.9em; }}
    main {{ max-width:1200px; margin:0 auto; padding:30px 20px; }}
    section {{ background:#fff; border-radius:8px; padding:24px; margin-bottom:24px;
               box-shadow:0 1px 4px rgba(0,0,0,0.08); }}
    h2 {{ border-bottom:2px solid #3498db; padding-bottom:6px; color:#2c3e50; margin-top:0; }}
    table {{ border-collapse:collapse; width:100%; font-size:0.88em; }}
    th {{ background:#3498db; color:#fff; padding:8px 10px; text-align:left; }}
    td {{ padding:7px 10px; border-bottom:1px solid #eee; }}
    tr:hover td {{ background:#f0f7ff; }}
    .badge {{ display:inline-block; padding:3px 10px; border-radius:12px;
              background:#e8f5e9; color:#2e7d32; font-size:0.8em; font-weight:600; }}
    .grid {{ display:grid; grid-template-columns:repeat(auto-fill,minmax(320px,1fr)); gap:16px; }}
    .card img {{ width:100%; border-radius:4px; }}
    .card p   {{ margin:6px 0 0; font-size:0.8em; color:#555; text-align:center; }}
    .stat-row {{ display:flex; gap:16px; flex-wrap:wrap; margin-bottom:12px; }}
    .stat {{ background:#f0f7ff; border-left:4px solid #3498db; padding:10px 16px;
             border-radius:4px; min-width:160px; }}
    .stat .val {{ font-size:1.5em; font-weight:700; color:#2c3e50; }}
    .stat .lbl {{ font-size:0.8em; color:#666; }}
  </style>
</head>
<body>
<header>
  <h1>🧬 scRNA-seq Pipeline Report</h1>
  <p>Generated on {date} &nbsp;|&nbsp; Pipeline: scRNA-seq v1.0</p>
</header>
<main>

<section>
  <h2>Dataset Overview</h2>
  <div class="stat-row">
    <div class="stat"><div class="val">{n_cells}</div><div class="lbl">Cells</div></div>
    <div class="stat"><div class="val">{n_genes}</div><div class="lbl">Genes</div></div>
    <div class="stat"><div class="val">{n_samples}</div><div class="lbl">Samples</div></div>
    <div class="stat"><div class="val">{n_clusters}</div><div class="lbl">Clusters</div></div>
  </div>
  {samples_table}
</section>

<section>
  <h2>Sample Sheet</h2>
  {sample_sheet_table}
</section>

<section>
  <h2>Clustering Summary</h2>
  {clustering_table}
</section>

<section>
  <h2>Top Marker Genes</h2>
  {markers_table}
</section>

<section>
  <h2>Figures</h2>
  <div class="grid">
    {figures_html}
  </div>
</section>

</main>
</body>
</html>
"""


def img_tag(path, caption):
    if not os.path.exists(path):
        return ""
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode()
    ext = os.path.splitext(path)[1].lstrip(".")
    return (f'<div class="card">'
            f'<img src="data:image/{ext};base64,{b64}" alt="{caption}">'
            f'<p>{caption}</p></div>')


def df_to_html(df, max_rows=20):
    return df.head(max_rows).to_html(index=True, classes="", border=0,
                                     float_format="{:.4f}".format)


def collect_figures(results_dir):
    figs = []
    for root, _, files in os.walk(results_dir):
        for fn in sorted(files):
            if fn.endswith(".png"):
                figs.append((os.path.join(root, fn),
                              os.path.relpath(os.path.join(root, fn), results_dir)))
    return figs


def main(args):
    from datetime import datetime

    adata   = ad.read_h5ad(args.annotated)
    samples = pd.read_csv(args.samples_tsv, sep="\t")
    clust   = pd.read_csv(args.clust_csv, index_col=0)
    markers = pd.read_csv(args.markers)

    n_cells    = f"{adata.n_obs:,}"
    n_genes    = f"{adata.n_vars:,}"
    n_samples  = str(len(samples))
    n_clusters = str(adata.obs["leiden"].nunique() if "leiden" in adata.obs.columns else "—")

    figs_html = "\n".join(
        img_tag(p, cap)
        for p, cap in collect_figures("results")
    )

    top_markers = (markers.groupby("cluster")
                   .apply(lambda x: x.nsmallest(3, "padj"))
                   .reset_index(drop=True)
                   [["cluster", "gene", "log2FC", "padj"]])

    html = TEMPLATE.format(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_cells=n_cells, n_genes=n_genes,
        n_samples=n_samples, n_clusters=n_clusters,
        samples_table="",
        sample_sheet_table=df_to_html(samples),
        clustering_table=df_to_html(clust),
        markers_table=df_to_html(top_markers),
        figures_html=figs_html,
    )

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Report saved: {args.out}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--annotated",   required=True)
    p.add_argument("--markers",     required=True)
    p.add_argument("--clust-csv",   required=True)
    p.add_argument("--samples-tsv", required=True)
    p.add_argument("--out",         default="results/report.html")
    main(p.parse_args())
