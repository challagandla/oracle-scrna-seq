#!/usr/bin/env python3
"""
Step 01 — Doublet Detection
============================
Identifies and removes doublets (two cells captured in one droplet).

Methods
-------
scrublet   (default) — pure Python, fast
scdblfinder            — R/Bioconductor, more accurate

Outputs
-------
results/<sample>/01_doublets/
  no_doublets.h5ad
  figures/
    01_doublet_score_umap.png   (if UMAP already computed)
    01_doublet_score_hist.png
"""

import argparse, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from utils.io_utils import load_adata, save_adata

sc.settings.verbosity = 1


def run_scrublet(adata, expected_rate=0.06):
    try:
        import scrublet as scr
    except ImportError:
        sys.exit("scrublet not installed. Run: pip install scrublet")

    import scipy.sparse as sp
    counts = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
    scrub = scr.Scrublet(counts, expected_doublet_rate=expected_rate)
    scores, predicted = scrub.scrub_doublets(
        min_counts=2, min_cells=3, n_prin_comps=30, verbose=False)
    adata.obs["doublet_score"]     = scores.astype(float)
    adata.obs["predicted_doublet"] = predicted
    rate = predicted.mean() * 100
    print(f"Scrublet: {predicted.sum()} doublets detected ({rate:.1f}%)")
    return scrub


def run_scdblfinder(adata, out_dir):
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
    except ImportError:
        sys.exit("rpy2 not installed. Run: conda install -c conda-forge rpy2")

    import pandas as pd
    tmp = os.path.join(out_dir, "_tmp_scdblfinder.h5ad")
    adata.write_h5ad(tmp)

    r_script = f"""
    suppressPackageStartupMessages({{
        library(scDblFinder); library(zellkonverter)
    }})
    sce  <- readH5AD("{tmp.replace(os.sep, '/')}")
    sce  <- scDblFinder(sce)
    df   <- as.data.frame(colData(sce)[, c("scDblFinder.score","scDblFinder.class")])
    write.csv(df, "{os.path.join(out_dir, '_scdblfinder.csv').replace(os.sep, '/')}")
    """
    ro.r(r_script)
    df = pd.read_csv(os.path.join(out_dir, "_scdblfinder.csv"), index_col=0)
    adata.obs["doublet_score"]     = df["scDblFinder.score"].values
    adata.obs["predicted_doublet"] = (df["scDblFinder.class"] == "doublet").values
    os.remove(tmp)
    n = int(adata.obs["predicted_doublet"].sum())
    print(f"scDblFinder: {n} doublets detected ({n/adata.n_obs*100:.1f}%)")


def plot_doublet_hist(adata, out_dir):
    fig, ax = plt.subplots(figsize=(5, 3))
    scores = adata.obs["doublet_score"]
    ax.hist(scores[~adata.obs["predicted_doublet"]], bins=50,
            alpha=0.6, color="steelblue", label="Singlet")
    ax.hist(scores[adata.obs["predicted_doublet"]],  bins=50,
            alpha=0.6, color="red",       label="Doublet")
    ax.set_xlabel("Doublet score"); ax.set_ylabel("Cells")
    ax.set_title("Doublet score distribution")
    ax.legend()
    plt.tight_layout()
    path = os.path.join(out_dir, "figures", "01_doublet_score_hist.png")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_adata(args.input)

    if args.method == "scrublet":
        run_scrublet(adata, expected_rate=args.expected_rate)
    else:
        run_scdblfinder(adata, args.out)

    plot_doublet_hist(adata, args.out)

    n_before = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    print(f"Cells after doublet removal: {adata.n_obs:,} "
          f"(removed {n_before - adata.n_obs:,})")

    save_adata(adata, os.path.join(args.out, "no_doublets.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",            default="results/doublets")
    p.add_argument("--method",         default="scrublet",
                   choices=["scrublet", "scdblfinder"])
    p.add_argument("--expected-rate",  type=float, default=0.06)
    main(p.parse_args())
