#!/usr/bin/env python3
"""Optional PAGA and RNA-velocity analysis with explicit alignment checks."""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from utils.io_utils import load_adata, save_adata

sc.settings.verbosity = 1


def run_paga(adata, cluster_key, fig_dir, seed=42):
    if cluster_key not in adata.obs:
        raise ValueError(f"PAGA cluster key '{cluster_key}' is missing")
    sc.tl.paga(adata, groups=cluster_key)
    sc.pl.paga(adata, show=False)
    plt.gcf().savefig(os.path.join(fig_dir, "01_paga_graph.png"), dpi=150, bbox_inches="tight")
    plt.close()
    sc.tl.umap(adata, init_pos="paga", random_state=seed)
    figure, axis = plt.subplots(figsize=(6, 5))
    sc.pl.umap(
        adata,
        color=cluster_key,
        ax=axis,
        show=False,
        frameon=False,
        title="PAGA-initialized UMAP",
        legend_loc="on data",
        legend_fontsize=7,
    )
    figure.savefig(os.path.join(fig_dir, "02_umap_paga_init.png"), dpi=150, bbox_inches="tight")
    plt.close(figure)


def _apply_barcode_map(loom, mapping_path):
    if not mapping_path:
        return loom
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
    required = {"analysis_barcode", "loom_barcode"}
    if set(mapping.columns) != required:
        raise ValueError(
            "barcode map must contain exactly analysis_barcode and loom_barcode columns"
        )
    for column in required:
        mapping[column] = mapping[column].str.strip()
    if (
        mapping.empty
        or mapping[list(required)].isna().any().any()
        or mapping[list(required)].eq("").any().any()
        or any(mapping[column].duplicated().any() for column in required)
    ):
        raise ValueError("barcode map must be complete and one-to-one")
    rename = mapping.set_index("loom_barcode")["analysis_barcode"]
    available = loom.obs_names.intersection(rename.index)
    if available.empty:
        raise ValueError("none of the mapped loom barcodes exists in the loom file")
    loom = loom[available].copy()
    loom.obs_names = rename.loc[loom.obs_names].to_numpy()
    return loom


def _validate_loom_identifiers(loom):
    loom.obs_names = loom.obs_names.astype(str)
    loom.var_names = loom.var_names.astype(str)
    if not loom.obs_names.is_unique:
        raise ValueError(
            "loom cell barcodes are duplicated; provide unique capture-level barcodes "
            "before velocity analysis"
        )
    if not loom.var_names.is_unique:
        raise ValueError(
            "loom gene identifiers are duplicated; resolve them against the count-matrix "
            "feature annotation before velocity analysis"
        )


def run_scvelo(
    adata,
    loom_path,
    cluster_key,
    fig_dir,
    *,
    barcode_map="",
    min_shared_cells=100,
    threads=1,
    seed=42,
):
    try:
        import scvelo as scv
    except ImportError as error:
        raise RuntimeError("scVelo is enabled but not installed") from error

    import numpy as np

    np.random.seed(seed)
    loom = sc.read_loom(loom_path)
    _validate_loom_identifiers(loom)
    loom = _apply_barcode_map(loom, barcode_map)

    shared_cells = adata.obs_names[adata.obs_names.isin(loom.obs_names)]
    shared_genes = adata.var_names[adata.var_names.isin(loom.var_names)]
    if len(shared_cells) < min_shared_cells:
        raise ValueError(
            f"only {len(shared_cells)} cells match the loom; expected at least "
            f"{min_shared_cells}. Supply trajectory.barcode_map_tsv when IDs differ."
        )
    if len(shared_genes) < 100:
        raise ValueError(f"only {len(shared_genes)} genes match the loom; check gene identifiers")

    velocity = adata[shared_cells, shared_genes].copy()
    aligned_loom = loom[shared_cells, shared_genes]
    for layer in ("spliced", "unspliced"):
        if layer not in aligned_loom.layers:
            raise ValueError(f"loom file is missing required '{layer}' layer")
        velocity.layers[layer] = aligned_loom.layers[layer].copy()
    if "ambiguous" in aligned_loom.layers:
        velocity.layers["ambiguous"] = aligned_loom.layers["ambiguous"].copy()

    scv.pp.filter_and_normalize(velocity, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(velocity, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(velocity, n_jobs=max(1, threads))
    scv.tl.velocity(velocity, mode="dynamical")
    scv.tl.velocity_graph(velocity, n_jobs=max(1, threads))
    scv.tl.latent_time(velocity)
    scv.tl.velocity_confidence(velocity)

    scv.pl.velocity_embedding_stream(
        velocity,
        basis="umap",
        color=cluster_key,
        show=False,
        save=False,
        figsize=(6, 5),
    )
    plt.savefig(os.path.join(fig_dir, "03_velocity_stream.png"), dpi=150, bbox_inches="tight")
    plt.close()
    scv.pl.scatter(velocity, color="latent_time", cmap="gnuplot", show=False, save=False)
    plt.savefig(os.path.join(fig_dir, "04_latent_time.png"), dpi=150, bbox_inches="tight")
    plt.close()
    scv.pl.scatter(
        velocity,
        color="velocity_confidence",
        cmap="RdYlGn",
        show=False,
        save=False,
    )
    plt.savefig(os.path.join(fig_dir, "05_velocity_confidence.png"), dpi=150, bbox_inches="tight")
    plt.close()
    return velocity


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    adata = load_adata(args.input)
    if args.paga:
        run_paga(adata, args.cluster_key, fig_dir, seed=args.seed)
    if args.loom:
        adata = run_scvelo(
            adata,
            args.loom,
            args.cluster_key,
            fig_dir,
            barcode_map=args.barcode_map,
            min_shared_cells=args.min_shared_cells,
            threads=args.threads,
            seed=args.seed,
        )
    save_adata(adata, os.path.join(args.out, "trajectory.h5ad"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--out", default="results/merged/09_trajectory")
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--paga", action="store_true")
    parser.add_argument("--loom", default="")
    parser.add_argument("--barcode-map", default="")
    parser.add_argument("--min-shared-cells", type=int, default=100)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    main(parser.parse_args())
