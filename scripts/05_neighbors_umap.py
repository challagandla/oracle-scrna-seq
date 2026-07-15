#!/usr/bin/env python3
"""
Step 05 — Neighbor Graph & UMAP
=================================
Builds the kNN graph used by both UMAP and Leiden, then embeds in 2D.

Outputs
-------
results/merged/05_embedding/
  embedded.h5ad
  figures/
    01_umap_qc_metrics.png
    02_umap_batch.png          (if batch_key provided)
    03_umap_sample.png
"""

import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import scanpy as sc

from utils.io_utils   import load_adata, save_adata
from utils.plot_utils import plot_umap_grid

sc.settings.verbosity = 1


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_adata(args.input)

    use_rep = args.use_rep
    if use_rep not in adata.obsm:
        raise ValueError(
            f"Requested representation '{use_rep}' is missing. Available: "
            f"{sorted(adata.obsm.keys())}. Check batch.method/embedding.use_rep."
        )

    n_pcs = min(args.n_pcs, adata.obsm[use_rep].shape[1])
    print(f"Building kNN graph: k={args.n_neighbors}, rep={use_rep}, pcs={args.n_pcs}")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors,
                    n_pcs=n_pcs, use_rep=use_rep,
                    metric=args.metric, random_state=args.seed)

    sc.tl.umap(adata, min_dist=args.min_dist, spread=args.spread,
               random_state=args.seed)
    print("UMAP computed")

    # QC overlays
    qc_cols = [c for c in ["total_counts", "n_genes_by_counts",
                            "pct_counts_mt", "pct_counts_ribo",
                            "doublet_score"] if c in adata.obs.columns]
    if qc_cols:
        plot_umap_grid(adata, qc_cols,
                       os.path.join(fig_dir, "01_umap_qc_metrics.png"),
                       title="UMAP — QC metrics")

    # Batch & sample overlays
    for col, fname in [
        (args.batch_key,  "02_umap_batch.png"),
        ("sample",        "03_umap_sample.png"),
        ("condition",     "04_umap_condition.png"),
    ]:
        if col and col in adata.obs.columns:
            plot_umap_grid(adata, [col],
                           os.path.join(fig_dir, fname))

    save_adata(adata, os.path.join(args.out, "embedded.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",          default="results/merged/05_embedding")
    p.add_argument("--n-neighbors",  type=int,   default=15)
    p.add_argument("--n-pcs",        type=int,   default=30)
    p.add_argument("--use-rep",      default="X_pca")
    p.add_argument("--metric",       default="euclidean")
    p.add_argument("--min-dist",     type=float, default=0.3)
    p.add_argument("--spread",       type=float, default=1.0)
    p.add_argument("--batch-key",    default="")
    p.add_argument("--seed",         type=int, default=42)
    main(p.parse_args())
