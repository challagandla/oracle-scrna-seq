#!/usr/bin/env python3
"""
Step 02 — Normalization
========================
Library-size normalization + log1p (default) or scran pooling.

Raw integer counts are always preserved in adata.layers["counts"].
adata.raw is frozen at the log-norm stage for downstream DE testing.

Outputs
-------
results/merged/02_normalize/
  normalized.h5ad
  figures/
    01_count_distribution_before.png
    02_count_distribution_after.png
    03_mean_variance_trend.png        (scran only)
"""

import argparse
import os
import sys
import tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.io_utils import load_adata, save_adata
from utils.validation import validate_raw_counts

sc.settings.verbosity = 1


def plot_count_dist(adata, title, path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    import scipy.sparse as sp
    data = adata.X
    cell_sums = (
        np.asarray(data.sum(axis=1)).ravel()
        if sp.issparse(data)
        else np.asarray(data).sum(axis=1)
    )
    axes[0].hist(cell_sums, bins=60, color="steelblue", edgecolor="none", alpha=0.8)
    axes[0].set_xlabel("Total counts per cell")
    axes[0].set_title("Library sizes")

    # Random sample of genes
    gene_means = np.asarray(data.mean(axis=0)).ravel()
    axes[1].hist(np.log1p(gene_means), bins=60, color="salmon", edgecolor="none", alpha=0.8)
    axes[1].set_xlabel("log(mean count + 1)")
    axes[1].set_title("Gene mean expression")

    fig.suptitle(title)
    plt.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def normalize_log(adata, target_sum=1e4):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.X = adata.X.astype(np.float32)
    adata.raw = adata  # freeze for DE testing
    print(f"Log-normalization complete (target_sum={target_sum:.0f})")


def normalize_scran(adata, seed=42):
    """Scran pooling-based size factors via rpy2."""
    try:
        import rpy2.robjects as ro
    except ImportError:
        sys.exit("rpy2 required for scran. Use --method log or install rpy2.")

    import scipy.sparse as sp

    adata.layers["counts"] = adata.X.copy()

    # Quick pre-clustering for pooling
    adata_tmp = adata.copy()
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)
    sc.pp.highly_variable_genes(adata_tmp, n_top_genes=2000)
    n_comps = min(20, adata_tmp.n_obs - 1, adata_tmp.n_vars - 1)
    if n_comps < 2:
        raise ValueError("scran requires at least 3 cells and 3 genes")
    sc.pp.pca(adata_tmp, n_comps=n_comps)
    sc.pp.neighbors(adata_tmp, n_neighbors=min(10, adata_tmp.n_obs - 1))
    sc.tl.leiden(adata_tmp, resolution=0.5, key_added="_precluster", random_state=seed)
    clusters = adata_tmp.obs["_precluster"].astype(int).values + 1

    # Let zellkonverter preserve the sparse matrix instead of materializing a
    # cells×genes dense array in Python/R memory.
    with tempfile.TemporaryDirectory(prefix="scran-") as tmp_dir:
        tmp_h5ad = os.path.join(tmp_dir, "counts.h5ad")
        adata.write_h5ad(tmp_h5ad)
        ro.globalenv["h5ad_path"] = ro.StrVector([tmp_h5ad])
        ro.globalenv["clusters_r"] = ro.IntVector(clusters.tolist())
        ro.globalenv["seed_r"] = ro.IntVector([seed])
        ro.r("""
            suppressPackageStartupMessages({
                library(scran)
                library(zellkonverter)
                library(SingleCellExperiment)
            })
            set.seed(seed_r[[1]])
            sce <- readH5AD(h5ad_path)
            if (!("counts" %in% assayNames(sce))) {
                counts(sce) <- assay(sce, "X")
            }
            size_factors <- as.numeric(calculateSumFactors(sce, clusters=clusters_r))
        """)
    sf = np.array(ro.r["size_factors"])
    if not np.isfinite(sf).all() or (sf <= 0).any():
        raise ValueError("scran returned non-finite or non-positive size factors")
    adata.obs["size_factor"] = sf

    if sp.issparse(adata.X):
        from scipy.sparse import diags
        adata.X = diags(1.0 / sf) @ adata.X
    else:
        adata.X = adata.X / sf[:, None]

    sc.pp.log1p(adata)
    adata.X = adata.X.astype(np.float32)
    adata.raw = adata
    print("Scran normalization complete")


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_adata(args.input)
    validate_raw_counts(adata, context=args.input)
    plot_count_dist(adata, "Before normalization",
                    os.path.join(fig_dir, "01_count_distribution_before.png"))

    if args.method == "log":
        normalize_log(adata, target_sum=args.target_sum)
    elif args.method == "scran":
        normalize_scran(adata, seed=args.seed)

    plot_count_dist(adata, "After normalization",
                    os.path.join(fig_dir, "02_count_distribution_after.png"))

    save_adata(adata, os.path.join(args.out, "normalized.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",         default="results/merged/02_normalize")
    p.add_argument("--method",      choices=["log", "scran"], default="log")
    p.add_argument("--target-sum",  type=float, default=1e4)
    p.add_argument("--seed",        type=int, default=42)
    main(p.parse_args())
