#!/usr/bin/env python3
"""
Step 06 — Leiden Clustering
=============================
Runs Leiden at multiple resolutions. Selects best by silhouette score.
Stores selected clustering in adata.obs["leiden"].

Outputs
-------
results/merged/06_cluster/
  clustered.h5ad
  clustering_summary.csv
  figures/
    01_umap_all_resolutions.png
    02_silhouette_curve.png
    03_umap_selected.png
    04_cluster_sizes.png
"""

import argparse, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.io_utils   import load_adata, save_adata
from utils.plot_utils import plot_umap_grid, plot_silhouette

sc.settings.verbosity = 1


def _silhouette(adata, cluster_key, use_rep):
    from sklearn.metrics import silhouette_score
    X   = adata.obsm[use_rep]
    lbl = adata.obs[cluster_key].astype(str).values
    if len(np.unique(lbl)) < 2:
        return np.nan
    n   = min(5000, X.shape[0])
    idx = np.random.default_rng(42).choice(X.shape[0], n, replace=False)
    try:
        return float(silhouette_score(X[idx], lbl[idx]))
    except Exception:
        return np.nan


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata    = load_adata(args.input)
    use_rep  = args.use_rep if args.use_rep in adata.obsm else "X_pca"
    rows     = []

    for res in args.resolutions:
        key = f"leiden_r{res:.2f}".replace(".", "_")
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_cl = adata.obs[key].nunique()
        sil  = _silhouette(adata, key, use_rep)
        rows.append({"resolution": res, "n_clusters": n_cl, "silhouette": sil})
        print(f"  res={res:.2f}  clusters={n_cl}  silhouette={sil:.4f}")

    df = pd.DataFrame(rows).set_index("resolution")
    df.to_csv(os.path.join(args.out, "clustering_summary.csv"))

    # Choose best resolution
    best_res = (df["silhouette"].idxmax()
                if not df["silhouette"].isna().all()
                else args.default_resolution)
    best_key = f"leiden_r{best_res:.2f}".replace(".", "_")
    adata.obs["leiden"] = adata.obs[best_key].copy()
    print(f"\nSelected resolution: {best_res}  "
          f"({int(df.loc[best_res, 'n_clusters'])} clusters)")

    # ── Plots ─────────────────────────────────────────────────
    # 1. All resolutions grid
    all_keys = [f"leiden_r{r:.2f}".replace(".", "_") for r in args.resolutions]
    plot_umap_grid(adata, all_keys,
                   os.path.join(fig_dir, "01_umap_all_resolutions.png"),
                   title="Leiden at all resolutions")

    # 2. Silhouette curve
    plot_silhouette(df, best_res,
                    os.path.join(fig_dir, "02_silhouette_curve.png"))

    # 3. Selected clustering
    plot_umap_grid(adata, ["leiden"],
                   os.path.join(fig_dir, "03_umap_selected.png"),
                   title=f"Selected clustering (res={best_res})")

    # 4. Cluster size bar
    counts = adata.obs["leiden"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(max(6, len(counts) * 0.5), 4))
    counts.plot.bar(ax=ax, color="steelblue")
    ax.set_xlabel("Cluster"); ax.set_ylabel("Cells")
    ax.set_title("Cluster sizes")
    plt.tight_layout()
    fig.savefig(os.path.join(fig_dir, "04_cluster_sizes.png"), dpi=150)
    plt.close(fig)

    save_adata(adata, os.path.join(args.out, "clustered.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",                default="results/merged/06_cluster")
    p.add_argument("--resolutions",        type=float, nargs="+",
                   default=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    p.add_argument("--default-resolution", type=float, default=0.6)
    p.add_argument("--use-rep",            default="X_pca")
    main(p.parse_args())
