#!/usr/bin/env python3
"""
Step 04 — Batch Correction
============================
Harmony (default), scVI, or ComBat.
If no --batch-key is provided the input is copied unchanged.

Outputs
-------
results/merged/04_batch/
  batch_corrected.h5ad
  figures/
    01_umap_before_correction.png   (colour = batch)
    02_umap_after_correction.png
"""

import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")

from utils.io_utils   import load_adata, save_adata
from utils.plot_utils import plot_umap_grid

sc.settings.verbosity = 1


def quick_umap(adata, use_rep, n_pcs, seed=42):
    """Compute a quick UMAP for visualisation."""
    dimensions = adata.obsm[use_rep].shape[1]
    sc.pp.neighbors(adata, use_rep=use_rep, n_pcs=min(n_pcs, dimensions))
    sc.tl.umap(adata, random_state=seed)


def harmony(adata, batch_key, seed):
    try:
        import harmonypy as hm
    except ImportError:
        sys.exit("harmonypy not installed. Run: pip install harmonypy")
    meta = adata.obs[[batch_key]]
    ho   = hm.run_harmony(adata.obsm["X_pca"], meta, batch_key,
                          max_iter_harmony=30, verbose=False,
                          random_state=seed)
    adata.obsm["X_pca_corrected"] = ho.Z_corr.T
    print("Harmony integration complete → X_pca_corrected")


def scvi_correct(adata, batch_key, n_latent, seed):
    try:
        import scvi
    except ImportError:
        sys.exit("scvi-tools not installed. Run: pip install scvi-tools")
    if "counts" not in adata.layers:
        sys.exit("scVI requires adata.layers['counts']. Ensure step 02 ran correctly.")
    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(max_epochs=400, early_stopping=True, plan_kwargs={"lr": 1e-3})
    adata.obsm["X_scVI"] = model.get_latent_representation()
    print(f"scVI training complete → X_scVI ({n_latent} dims)")


def combat(adata, batch_key, seed):
    # Correct only HVGs. Full-matrix ComBat is needlessly dense for scRNA data
    # and downstream neighbors use the corrected PCA representation, not X.
    hvg_mask = adata.var.get("highly_variable", np.ones(adata.n_vars, dtype=bool))
    corrected_hvg = adata[:, np.asarray(hvg_mask, dtype=bool)].copy()
    sc.pp.combat(corrected_hvg, key=batch_key)
    sc.pp.scale(corrected_hvg, max_value=10)
    n_comps = min(50, corrected_hvg.n_obs - 1, corrected_hvg.n_vars - 1)
    if n_comps < 2:
        raise ValueError("ComBat PCA requires at least 3 cells and 3 selected genes")
    sc.tl.pca(
        corrected_hvg,
        n_comps=n_comps,
        svd_solver="arpack",
        random_state=seed,
    )
    adata.obsm["X_pca_corrected"] = corrected_hvg.obsm["X_pca"].copy()
    print("ComBat correction applied → X_pca_corrected")


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_adata(args.input)

    if args.method == "none" or not args.batch_key:
        reason = "method=none" if args.method == "none" else "no batch key provided"
        print(f"Skipping batch correction ({reason}).")
        save_adata(adata, os.path.join(args.out, "batch_corrected.h5ad"))
        return

    if args.batch_key not in adata.obs:
        raise ValueError(f"Batch key '{args.batch_key}' is absent from AnnData metadata")
    n_batches = adata.obs[args.batch_key].nunique()
    if n_batches < 2:
        raise ValueError("Batch correction requires at least two batch levels")
    print(f"Batch key: '{args.batch_key}'  |  {n_batches} batches")

    # Pre-correction UMAP
    if "X_umap" not in adata.obsm:
        quick_umap(adata, use_rep="X_pca", n_pcs=30, seed=args.seed)
    plot_umap_grid(adata, [args.batch_key],
                   os.path.join(fig_dir, "01_umap_before_correction.png"),
                   title="Before batch correction")

    if args.method == "harmony":
        harmony(adata, args.batch_key, args.seed)
        use_rep = "X_pca_corrected"
    elif args.method == "scvi":
        scvi_correct(adata, args.batch_key, args.n_latent, args.seed)
        use_rep = "X_scVI"
    elif args.method == "combat":
        combat(adata, args.batch_key, args.seed)
        use_rep = "X_pca_corrected"

    # Post-correction UMAP
    quick_umap(adata, use_rep=use_rep, n_pcs=30, seed=args.seed)
    plot_umap_grid(adata, [args.batch_key],
                   os.path.join(fig_dir, "02_umap_after_correction.png"),
                   title=f"After batch correction ({args.method})")

    save_adata(adata, os.path.join(args.out, "batch_corrected.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",        default="results/merged/04_batch")
    p.add_argument("--method",     choices=["none", "harmony", "scvi", "combat"],
                   default="harmony")
    p.add_argument("--batch-key",  default="")
    p.add_argument("--n-latent",   type=int, default=30)
    p.add_argument("--seed",       type=int, default=42)
    main(p.parse_args())
