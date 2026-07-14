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

import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from utils.io_utils import load_adata, save_adata
from utils.validation import validate_raw_counts

sc.settings.verbosity = 1


def run_scrublet(adata, expected_rate=0.06, seed=42):
    try:
        import scrublet as scr
    except ImportError:
        sys.exit("scrublet not installed. Run: pip install scrublet")

    if adata.n_obs < 30:
        raise ValueError("Scrublet requires at least 30 cells in each capture/library")
    n_components = min(30, adata.n_obs - 1, adata.n_vars - 1)
    counts = adata.X.copy()
    scrub = scr.Scrublet(
        counts,
        expected_doublet_rate=expected_rate,
        random_state=seed,
    )
    scores, predicted = scrub.scrub_doublets(
        min_counts=2, min_cells=3, n_prin_comps=n_components, verbose=False)
    adata.obs["doublet_score"]     = scores.astype(float)
    adata.obs["predicted_doublet"] = predicted
    rate = predicted.mean() * 100
    print(f"Scrublet: {predicted.sum()} doublets detected ({rate:.1f}%)")
    return scrub


def run_scdblfinder(adata, out_dir, expected_rate=0.06, seed=42):
    try:
        import rpy2.robjects as ro
    except ImportError:
        sys.exit("rpy2 not installed. Run: conda install -c conda-forge rpy2")

    import pandas as pd
    tmp = os.path.join(out_dir, "_tmp_scdblfinder.h5ad")
    adata.write_h5ad(tmp)

    csv_path = os.path.join(out_dir, "_scdblfinder.csv")
    ro.globalenv["input_h5ad"] = ro.StrVector([tmp.replace(os.sep, "/")])
    ro.globalenv["output_csv"] = ro.StrVector([csv_path.replace(os.sep, "/")])
    ro.globalenv["seed_r"] = ro.IntVector([seed])
    ro.globalenv["expected_rate_r"] = ro.FloatVector([expected_rate])
    r_script = """
    suppressPackageStartupMessages({
        library(scDblFinder); library(zellkonverter); library(SingleCellExperiment)
    })
    set.seed(seed_r[[1]])
    sce  <- readH5AD(input_h5ad)
    if (!("counts" %in% assayNames(sce))) counts(sce) <- assay(sce, "X")
    sce  <- scDblFinder(sce, dbr=expected_rate_r[[1]])
    df   <- as.data.frame(colData(sce)[, c("scDblFinder.score","scDblFinder.class")])
    write.csv(df, output_csv)
    """
    try:
        ro.r(r_script)
        df = pd.read_csv(csv_path, index_col=0).reindex(adata.obs_names)
        if df.isna().any().any():
            raise ValueError("scDblFinder output barcodes do not match the input AnnData")
        adata.obs["doublet_score"] = df["scDblFinder.score"].to_numpy()
        adata.obs["predicted_doublet"] = (
            df["scDblFinder.class"] == "doublet"
        ).to_numpy()
    finally:
        for path in (tmp, csv_path):
            if os.path.exists(path):
                os.remove(path)
    n = int(adata.obs["predicted_doublet"].sum())
    print(f"scDblFinder: {n} doublets detected ({n/adata.n_obs*100:.1f}%)")


def plot_doublet_hist(adata, out_dir):
    fig, ax = plt.subplots(figsize=(5, 3))
    scores = adata.obs["doublet_score"]
    ax.hist(scores[~adata.obs["predicted_doublet"]], bins=50,
            alpha=0.6, color="steelblue", label="Singlet")
    ax.hist(scores[adata.obs["predicted_doublet"]],  bins=50,
            alpha=0.6, color="red",       label="Doublet")
    ax.set_xlabel("Doublet score")
    ax.set_ylabel("Cells")
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
    validate_raw_counts(adata, context=args.input)

    if args.method == "scrublet":
        run_scrublet(adata, expected_rate=args.expected_rate, seed=args.seed)
    else:
        run_scdblfinder(
            adata,
            args.out,
            expected_rate=args.expected_rate,
            seed=args.seed,
        )

    plot_doublet_hist(adata, args.out)

    n_before = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    if adata.n_obs == 0:
        raise ValueError("Doublet detection removed every cell; review the capture and threshold")
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
    p.add_argument("--seed",           type=int, default=42)
    main(p.parse_args())
