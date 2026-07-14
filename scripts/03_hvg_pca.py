#!/usr/bin/env python3
"""
Step 03 — HVG Selection & PCA
==============================
Selects highly variable genes (HVGs) and computes PCA.

Outputs
-------
results/merged/03_hvg_pca/
  hvg_pca.h5ad
  figures/
    01_hvg_dispersion.png
    02_pca_elbow.png
    03_pca_scatter_qc.png      PC1 vs PC2 coloured by QC metrics
    04_pca_loadings.png        top contributing genes per PC
"""

import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.io_utils   import load_adata, save_adata
from utils.plot_utils import plot_elbow

sc.settings.verbosity = 1


def select_hvg(adata, n_top, flavor, batch_key=None):
    """seurat_v3 expects raw counts; swap in layers['counts'] if available."""
    n_top = min(int(n_top), adata.n_vars)
    if flavor == "seurat_v3" and "counts" in adata.layers:
        orig_X = adata.X
        adata.X = adata.layers["counts"].copy()
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=flavor,
                                        batch_key=batch_key or None, subset=False)
        finally:
            adata.X = orig_X
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=flavor,
                                    batch_key=batch_key or None, subset=False)
    n = int(adata.var["highly_variable"].sum())
    print(f"HVG selected: {n} / {adata.n_vars} genes  (flavor={flavor})")


def plot_hvg_dispersion(adata, path):
    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    fig.set_size_inches(6, 4)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_pca(adata, n_comps, seed=42):
    hvg_adata = adata[:, adata.var["highly_variable"]].copy()
    n_comps = min(int(n_comps), hvg_adata.n_obs - 1, hvg_adata.n_vars - 1)
    if n_comps < 2:
        raise ValueError("PCA requires at least 3 cells and 3 selected HVGs")
    sc.pp.scale(hvg_adata, max_value=10)
    sc.tl.pca(hvg_adata, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    # Copy results back to full adata
    adata.obsm["X_pca"] = hvg_adata.obsm["X_pca"]
    adata.varm["PCs"]   = np.zeros((adata.n_vars, n_comps))
    hvg_idx = np.where(adata.var["highly_variable"])[0]
    adata.varm["PCs"][hvg_idx] = hvg_adata.varm["PCs"]
    adata.uns["pca"]    = hvg_adata.uns["pca"]
    print(f"PCA: {n_comps} components computed")
    return n_comps


def plot_pca_qc(adata, path):
    qc_cols = [c for c in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
               if c in adata.obs.columns]
    n = len(qc_cols)
    if n == 0:
        return
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), squeeze=False)
    for ax, c in zip(axes[0], qc_cols, strict=True):
        sc.pl.pca(adata, color=c, ax=ax, show=False, title=c)
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_pca_loadings(adata, n_top, path):
    """Bar plot of top genes contributing to PC1 and PC2."""
    PCs = adata.varm["PCs"]
    genes = adata.var_names.to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for pc_idx, ax in enumerate(axes):
        loadings = PCs[:, pc_idx]
        top_idx  = np.argsort(np.abs(loadings))[::-1][:n_top]
        ax.barh(genes[top_idx][::-1], loadings[top_idx][::-1],
                color=["red" if v > 0 else "steelblue" for v in loadings[top_idx][::-1]])
        ax.set_title(f"PC{pc_idx+1} top {n_top} loadings")
        ax.axvline(0, color="black", lw=0.8)
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_adata(args.input)

    select_hvg(adata, n_top=args.n_hvg, flavor=args.hvg_flavor,
               batch_key=args.batch_key or None)
    plot_hvg_dispersion(adata, os.path.join(fig_dir, "01_hvg_dispersion.png"))

    n_comps = run_pca(adata, n_comps=args.n_pcs, seed=args.seed)
    plot_elbow(adata, n_pcs_use=min(args.n_pcs_use, n_comps),
               path=os.path.join(fig_dir, "02_pca_elbow.png"))
    plot_pca_qc(adata, os.path.join(fig_dir, "03_pca_scatter_qc.png"))
    plot_pca_loadings(adata, n_top=15,
                      path=os.path.join(fig_dir, "04_pca_loadings.png"))

    save_adata(adata, os.path.join(args.out, "hvg_pca.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",          default="results/merged/03_hvg_pca")
    p.add_argument("--n-hvg",        type=int, default=2000)
    p.add_argument("--hvg-flavor",   default="seurat_v3",
                   choices=["seurat_v3", "seurat", "cell_ranger"])
    p.add_argument("--n-pcs",        type=int, default=50)
    p.add_argument("--n-pcs-use",    type=int, default=30,
                   help="PCs to use downstream (marked on elbow plot)")
    p.add_argument("--batch-key",    default="")
    p.add_argument("--seed",         type=int, default=42)
    main(p.parse_args())
