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

import argparse, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.io_utils import load_adata, save_adata

sc.settings.verbosity = 1


def plot_count_dist(adata, title, path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    import scipy.sparse as sp
    data = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
    cell_sums = data.sum(axis=1)
    axes[0].hist(cell_sums, bins=60, color="steelblue", edgecolor="none", alpha=0.8)
    axes[0].set_xlabel("Total counts per cell"); axes[0].set_title("Library sizes")

    # Random sample of genes
    gene_means = data.mean(axis=0)
    axes[1].hist(np.log1p(gene_means), bins=60, color="salmon", edgecolor="none", alpha=0.8)
    axes[1].set_xlabel("log(mean count + 1)"); axes[1].set_title("Gene mean expression")

    fig.suptitle(title)
    plt.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def normalize_log(adata, target_sum=1e4):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata  # freeze for DE testing
    print(f"Log-normalization complete (target_sum={target_sum:.0f})")


def normalize_scran(adata):
    """Scran pooling-based size factors via rpy2."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri
        numpy2ri.activate()
    except ImportError:
        sys.exit("rpy2 required for scran. Use --method log or install rpy2.")

    import scipy.sparse as sp

    adata.layers["counts"] = adata.X.copy()

    # Quick pre-clustering for pooling
    adata_tmp = adata.copy()
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)
    sc.pp.highly_variable_genes(adata_tmp, n_top_genes=2000)
    sc.pp.pca(adata_tmp, n_comps=20)
    sc.pp.neighbors(adata_tmp, n_neighbors=10)
    sc.tl.leiden(adata_tmp, resolution=0.5, key_added="_precluster")
    clusters = adata_tmp.obs["_precluster"].astype(int).values + 1

    mat = adata.X.toarray() if sp.issparse(adata.X) else adata.X.astype(float)
    ro.r.assign("counts_r",   mat.T)
    ro.r.assign("clusters_r", clusters)
    ro.r("""
        suppressPackageStartupMessages(library(scran))
        sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts_r))
        sf  <- calculateSumFactors(sce, clusters=clusters_r)
        size_factors <- as.numeric(sf)
    """)
    sf = np.array(ro.r["size_factors"])
    adata.obs["size_factor"] = sf

    if sp.issparse(adata.X):
        from scipy.sparse import diags
        adata.X = diags(1.0 / sf) @ adata.X
    else:
        adata.X = adata.X / sf[:, None]

    sc.pp.log1p(adata)
    adata.raw = adata
    print("Scran normalization complete")


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_adata(args.input)
    plot_count_dist(adata, "Before normalization",
                    os.path.join(fig_dir, "01_count_distribution_before.png"))

    if args.method == "log":
        normalize_log(adata, target_sum=args.target_sum)
    elif args.method == "scran":
        normalize_scran(adata)

    plot_count_dist(adata, "After normalization",
                    os.path.join(fig_dir, "02_count_distribution_after.png"))

    save_adata(adata, os.path.join(args.out, "normalized.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",         default="results/merged/02_normalize")
    p.add_argument("--method",      choices=["log", "scran"], default="log")
    p.add_argument("--target-sum",  type=float, default=1e4)
    main(p.parse_args())
