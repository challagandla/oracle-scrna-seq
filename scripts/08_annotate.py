#!/usr/bin/env python3
"""
Step 08 — Cell Type Annotation
================================
Three complementary approaches (use in combination):
  1. Gene-score-based   — sc.tl.score_genes with curated marker lists
  2. CellTypist         — ML classifier (immune cells / PBMC atlas)
  3. Manual             — load cluster→cell_type TSV

Outputs
-------
results/merged/08_annotate/
  annotated.h5ad
  figures/
    01_umap_leiden.png
    02_umap_marker_scores.png
    03_umap_predicted_celltype.png
    04_umap_celltypist.png         (if --celltypist)
    05_umap_manual.png             (if --manual-tsv)
    06_composition_barplot.png
"""

import argparse, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.io_utils   import load_adata, save_adata
from utils.plot_utils import plot_umap_grid

sc.settings.verbosity = 1


# ── Built-in marker dictionary (human PBMC) ──────────────────────────────────
BUILTIN_MARKERS = {
    "T cell":          ["CD3D", "CD3E", "CD3G"],
    "CD4 T cell":      ["CD4", "IL7R", "CCR7", "SELL"],
    "CD8 T cell":      ["CD8A", "CD8B", "GZMB", "GZMK"],
    "Treg":            ["FOXP3", "IL2RA", "CTLA4"],
    "NK cell":         ["GNLY", "NKG7", "NCAM1", "KLRD1"],
    "B cell":          ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "Plasma cell":     ["MZB1", "SDC1", "JCHAIN", "IGHG1"],
    "Classical Mono":  ["CD14", "LYZ", "S100A8", "S100A9"],
    "Non-classical Mono": ["FCGR3A", "MS4A7", "CX3CR1"],
    "mDC":             ["FCER1A", "CST3", "CLEC9A"],
    "pDC":             ["LILRA4", "IL3RA", "CLEC4C"],
    "Platelet":        ["PPBP", "PF4", "GP1BA"],
    "Erythrocyte":     ["HBB", "HBA1", "HBA2"],
}

MOUSE_MARKERS = {
    "T cell":          ["Cd3d", "Cd3e", "Cd3g"],
    "CD4 T cell":      ["Cd4", "Il7r", "Ccr7"],
    "CD8 T cell":      ["Cd8a", "Cd8b1", "Gzmb"],
    "NK cell":         ["Klrb1c", "Nkg7", "Ncr1"],
    "B cell":          ["Cd19", "Ms4a1", "Cd79a"],
    "Macrophage":      ["Adgre1", "Cd68", "Csf1r"],
    "Microglia":       ["P2ry12", "Tmem119", "Cx3cr1"],
    "Neuron":          ["Rbfox3", "Map2", "Tubb3"],
    "Astrocyte":       ["Gfap", "Aqp4", "Aldh1l1"],
    "Oligodendrocyte": ["Mbp", "Mog", "Plp1"],
}


def score_markers(adata, marker_dict, cluster_key):
    score_keys = []
    for ct, genes in marker_dict.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        key = f"score_{ct.replace(' ', '_').replace('/', '_')}"
        sc.tl.score_genes(adata, gene_list=present,
                          score_name=key,
                          use_raw=(adata.raw is not None))
        score_keys.append((ct, key))

    if not score_keys:
        print("No marker genes found in dataset")
        return

    # Assign best cell type per cluster
    score_cols = [k for _, k in score_keys]
    cluster_mean = adata.obs.groupby(cluster_key)[score_cols].mean()
    cluster_mean.columns = [ct for ct, _ in score_keys]
    best = cluster_mean.idxmax(axis=1)
    adata.obs["predicted_celltype"] = (
        adata.obs[cluster_key].map(best).astype("category")
    )
    print("\nCluster → predicted cell type:")
    for cl, ct in best.items():
        print(f"  Cluster {cl:>4} → {ct}")
    return score_keys


def run_celltypist(adata, model_name):
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        sys.exit("Install CellTypist:  pip install celltypist")

    models.download_models(model=model_name, force_update=False)
    model = models.Model.load(model=model_name)
    pred  = celltypist.annotate(adata, model=model, majority_voting=True)
    adata.obs["celltypist_type"]     = pred.predicted_labels["predicted_labels"]
    adata.obs["celltypist_conf"]     = pred.predicted_labels["conf_score"]
    adata.obs["celltypist_majority"] = pred.predicted_labels["majority_voting"]
    print("CellTypist annotation complete")


def load_manual(adata, tsv_path, cluster_key):
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).set_index("cluster")
    if "cell_type" not in df.columns:
        sys.exit("Manual TSV must have columns: cluster, cell_type")
    adata.obs["cell_type_manual"] = (
        adata.obs[cluster_key].map(df["cell_type"].to_dict()).astype("category")
    )
    print("Manual annotations loaded")


def composition_barplot(adata, group_key, split_key, path):
    """Stacked bar plot: % of each cell type per sample/batch."""
    if group_key not in adata.obs.columns or split_key not in adata.obs.columns:
        return
    ct = pd.crosstab(adata.obs[split_key], adata.obs[group_key], normalize="index")
    fig, ax = plt.subplots(figsize=(max(8, len(ct.index)), 5))
    ct.plot.bar(stacked=True, ax=ax, colormap="tab20", legend=True, width=0.8)
    ax.set_ylabel("Proportion"); ax.set_xlabel(split_key)
    ax.set_title(f"Cell type composition per {split_key}")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)
    plt.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {os.path.basename(path)}")


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata       = load_adata(args.input)
    cluster_key = args.cluster_key

    # ── UMAP with cluster labels ──────────────────────────────
    plot_umap_grid(adata, [cluster_key],
                   os.path.join(fig_dir, "01_umap_leiden.png"))

    # ── Marker scoring ────────────────────────────────────────
    if args.markers_json:
        with open(args.markers_json) as f:
            marker_dict = json.load(f)
    else:
        marker_dict = BUILTIN_MARKERS

    score_keys = score_markers(adata, marker_dict, cluster_key)
    if score_keys:
        score_cols = [k for _, k in score_keys[:6]]  # max 6 for grid
        plot_umap_grid(adata, score_cols,
                       os.path.join(fig_dir, "02_umap_marker_scores.png"),
                       title="Gene module scores")
        plot_umap_grid(adata, ["predicted_celltype"],
                       os.path.join(fig_dir, "03_umap_predicted_celltype.png"))

    # ── CellTypist ────────────────────────────────────────────
    if args.celltypist:
        run_celltypist(adata, model_name=args.celltypist_model)
        plot_umap_grid(adata, ["celltypist_type", "celltypist_majority"],
                       os.path.join(fig_dir, "04_umap_celltypist.png"))

    # ── Manual ────────────────────────────────────────────────
    if args.manual_tsv:
        load_manual(adata, args.manual_tsv, cluster_key)
        plot_umap_grid(adata, ["cell_type_manual"],
                       os.path.join(fig_dir, "05_umap_manual.png"))

    # ── Composition plot ──────────────────────────────────────
    annot_key = ("cell_type_manual" if args.manual_tsv else
                 "celltypist_majority" if args.celltypist else
                 "predicted_celltype")
    if "sample" in adata.obs.columns:
        composition_barplot(adata, annot_key, "sample",
                            os.path.join(fig_dir, "06_composition_barplot.png"))

    save_adata(adata, os.path.join(args.out, "annotated.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",              default="results/merged/08_annotate")
    p.add_argument("--cluster-key",      default="leiden")
    p.add_argument("--markers-json",     default="")
    p.add_argument("--celltypist",       action="store_true")
    p.add_argument("--celltypist-model", default="Immune_All_Low.pkl")
    p.add_argument("--manual-tsv",       default="")
    main(p.parse_args())
