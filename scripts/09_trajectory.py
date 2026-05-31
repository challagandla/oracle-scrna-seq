#!/usr/bin/env python3
"""
Step 09 — Trajectory Analysis (optional)
==========================================
PAGA   — cluster-level connectivity graph (no splicing needed)
scVelo — RNA velocity (requires .loom with spliced/unspliced counts)

Outputs
-------
results/merged/09_trajectory/
  trajectory.h5ad
  figures/
    01_paga_graph.png
    02_umap_paga_init.png
    03_velocity_stream.png   (scVelo)
    04_latent_time.png       (scVelo)
    05_velocity_confidence.png
"""

import argparse, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.io_utils import load_adata, save_adata

sc.settings.verbosity = 1


def run_paga(adata, cluster_key, fig_dir):
    sc.tl.paga(adata, groups=cluster_key)

    # PAGA graph
    sc.pl.paga(adata, show=False)
    fig = plt.gcf()
    fig.savefig(os.path.join(fig_dir, "01_paga_graph.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)

    # PAGA-initialized UMAP
    sc.tl.umap(adata, init_pos="paga")
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata, color=cluster_key, ax=ax2, show=False,
               frameon=False, title="PAGA-initialized UMAP",
               legend_loc="on data", legend_fontsize=7)
    fig2.savefig(os.path.join(fig_dir, "02_umap_paga_init.png"),
                 dpi=150, bbox_inches="tight")
    plt.close(fig2)

    # PAGA abstracted graph with threshold
    sc.pl.paga(adata, threshold=0.05, show=False)
    fig3 = plt.gcf()
    fig3.savefig(os.path.join(fig_dir, "01b_paga_graph_threshold.png"),
                 dpi=150, bbox_inches="tight")
    plt.close(fig3)
    print("PAGA complete")


def run_scvelo(adata, loom_path, cluster_key, fig_dir):
    try:
        import scvelo as scv
    except ImportError:
        sys.exit("scvelo not installed. Run: pip install scvelo")

    ldata = sc.read_loom(loom_path)
    ldata.var_names_make_unique()

    # Match barcodes
    shared = adata.obs_names.intersection(ldata.obs_names)
    if len(shared) < 100:
        sys.exit(f"Only {len(shared)} shared barcodes. Check loom file.")
    print(f"Shared barcodes: {len(shared)}")

    adata_v = adata[shared].copy()
    adata_v.layers["spliced"]   = ldata[shared].layers["spliced"]
    adata_v.layers["unspliced"] = ldata[shared].layers["unspliced"]
    adata_v.layers["ambiguous"] = ldata[shared].layers.get("ambiguous", None)

    scv.pp.filter_and_normalize(adata_v, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata_v, n_pcs=30, n_neighbors=30)

    # Dynamical model
    scv.tl.recover_dynamics(adata_v, n_jobs=-1)
    scv.tl.velocity(adata_v, mode="dynamical")
    scv.tl.velocity_graph(adata_v)
    scv.tl.latent_time(adata_v)
    scv.tl.velocity_confidence(adata_v)

    # Plots
    scv.pl.velocity_embedding_stream(
        adata_v, basis="umap", color=cluster_key,
        show=False, save=False, figsize=(6, 5))
    plt.savefig(os.path.join(fig_dir, "03_velocity_stream.png"),
                dpi=150, bbox_inches="tight")
    plt.close()

    scv.pl.scatter(adata_v, color="latent_time",
                   cmap="gnuplot", show=False, save=False)
    plt.savefig(os.path.join(fig_dir, "04_latent_time.png"),
                dpi=150, bbox_inches="tight")
    plt.close()

    scv.pl.scatter(adata_v, color="velocity_confidence",
                   cmap="RdYlGn", show=False, save=False)
    plt.savefig(os.path.join(fig_dir, "05_velocity_confidence.png"),
                dpi=150, bbox_inches="tight")
    plt.close()

    # Save velocity subset
    save_adata(adata_v, os.path.join(os.path.dirname(fig_dir), "velocity.h5ad"))
    print("scVelo complete")
    return adata_v


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_adata(args.input)

    if args.paga:
        print("\n── PAGA ──")
        run_paga(adata, args.cluster_key, fig_dir)

    if args.loom:
        print("\n── scVelo ──")
        run_scvelo(adata, args.loom, args.cluster_key, fig_dir)

    save_adata(adata, os.path.join(args.out, "trajectory.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",          default="results/merged/09_trajectory")
    p.add_argument("--cluster-key",  default="leiden")
    p.add_argument("--paga",         action="store_true")
    p.add_argument("--loom",         default="")
    main(p.parse_args())
