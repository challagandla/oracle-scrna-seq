#!/usr/bin/env python3
"""Conservative, opt-in cell-type annotation with explicit uncertainty."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from utils.io_utils import load_adata, save_adata
from utils.plot_utils import plot_umap_grid

sc.settings.verbosity = 1

MARKER_SETS = {
    "human_pbmc": {
        "T cell": ["CD3D", "CD3E", "CD3G"],
        "CD4 T cell": ["CD4", "IL7R", "CCR7", "SELL"],
        "CD8 T cell": ["CD8A", "CD8B", "GZMK", "NKG7"],
        "NK cell": ["GNLY", "NKG7", "NCAM1", "KLRD1"],
        "B cell": ["MS4A1", "CD79A", "CD79B", "CD37"],
        "Plasma cell": ["MZB1", "SDC1", "JCHAIN", "IGHG1"],
        "Classical monocyte": ["CD14", "LYZ", "S100A8", "S100A9"],
        "Non-classical monocyte": ["FCGR3A", "MS4A7", "LST1", "CX3CR1"],
        "Dendritic cell": ["FCER1A", "CST3", "CLEC10A"],
        "Platelet": ["PPBP", "PF4", "GP1BA"],
    },
    "mouse_brain": {
        "T cell": ["Cd3d", "Cd3e", "Cd3g"],
        "B cell": ["Cd19", "Ms4a1", "Cd79a"],
        "Macrophage": ["Adgre1", "Cd68", "Csf1r"],
        "Microglia": ["P2ry12", "Tmem119", "Cx3cr1"],
        "Neuron": ["Rbfox3", "Map2", "Tubb3"],
        "Astrocyte": ["Gfap", "Aqp4", "Aldh1l1"],
        "Oligodendrocyte": ["Mbp", "Mog", "Plp1"],
    },
}


def load_marker_dictionary(path):
    """Load and validate a custom marker dictionary."""
    with open(path, encoding="utf-8") as handle:
        marker_dictionary = json.load(handle)
    if not isinstance(marker_dictionary, dict) or not marker_dictionary:
        raise ValueError("custom marker JSON must be a non-empty object")
    validated = {}
    for label, genes in marker_dictionary.items():
        clean_label = str(label).strip()
        if not clean_label or not isinstance(genes, list):
            raise ValueError("each marker label must be non-blank and map to a gene list")
        clean_genes = list(
            dict.fromkeys(str(gene).strip() for gene in genes if str(gene).strip())
        )
        if len(clean_genes) < 2:
            raise ValueError(f"marker label '{clean_label}' needs at least two unique genes")
        if clean_label in validated:
            raise ValueError(f"marker label '{clean_label}' is duplicated after trimming")
        validated[clean_label] = clean_genes
    return validated


def score_markers(adata, marker_dict, cluster_key, min_coverage, min_margin, seed=42):
    """Assign provisional labels only when marker coverage and margin are adequate."""
    score_keys = []
    coverage_rows = []
    for cell_type, genes in marker_dict.items():
        unique_genes = list(dict.fromkeys(genes))
        present = [gene for gene in unique_genes if gene in adata.var_names]
        required = max(2, math.ceil(len(unique_genes) * min_coverage))
        coverage_rows.append(
            {
                "cell_type": cell_type,
                "markers_present": len(present),
                "markers_total": len(unique_genes),
                "eligible": len(present) >= required,
            }
        )
        if len(present) < required:
            continue
        key = f"marker_score_{len(score_keys)}"
        sc.tl.score_genes(
            adata,
            gene_list=present,
            score_name=key,
            use_raw=adata.raw is not None,
            random_state=seed,
        )
        score_keys.append((cell_type, key))

    if not score_keys:
        raise ValueError(
            "No cell type met the marker-coverage threshold. Check species, gene "
            "identifiers, tissue, and annotation.marker_set/markers_json."
        )

    score_columns = [key for _, key in score_keys]
    cluster_means = adata.obs.groupby(cluster_key, observed=True)[score_columns].mean()
    cluster_means.columns = [cell_type for cell_type, _ in score_keys]
    summary_rows = []
    labels = {}
    for cluster, scores in cluster_means.iterrows():
        ordered = scores.sort_values(ascending=False)
        best_type = str(ordered.index[0])
        best_score = float(ordered.iloc[0])
        second_score = float(ordered.iloc[1]) if len(ordered) > 1 else 0.0
        margin = best_score - second_score
        label = best_type if best_score > 0 and margin >= min_margin else "Unknown"
        labels[cluster] = label
        summary_rows.append(
            {
                "cluster": str(cluster),
                "marker_label": label,
                "best_marker_type": best_type,
                "best_score": best_score,
                "score_margin": margin,
            }
        )

    adata.obs["marker_cell_type"] = adata.obs[cluster_key].map(labels).astype("category")
    return score_keys, pd.DataFrame(summary_rows), pd.DataFrame(coverage_rows)


def prepare_celltypist_input(adata):
    """Build CellTypist's required log1p, 10,000-count input from raw counts."""
    if "counts" not in adata.layers:
        raise ValueError(
            "CellTypist requires layers['counts']; rerun normalization from raw counts"
        )
    annotation = ad.AnnData(
        X=adata.layers["counts"].copy(),
        obs=pd.DataFrame(index=adata.obs_names.copy()),
        var=pd.DataFrame(index=adata.var_names.copy()),
    )
    sc.pp.normalize_total(annotation, target_sum=10_000)
    sc.pp.log1p(annotation)
    return annotation


def run_celltypist(adata, model_name, min_confidence):
    """Run CellTypist and preserve an Unknown label below the confidence floor."""
    try:
        import celltypist
        from celltypist import models
    except ImportError as error:
        raise RuntimeError("CellTypist is enabled but not installed") from error

    models.download_models(model=model_name, force_update=False)
    model = models.Model.load(model=model_name)
    annotation = prepare_celltypist_input(adata)
    prediction = celltypist.annotate(annotation, model=model, majority_voting=True)
    labels = prediction.predicted_labels.reindex(adata.obs_names)
    if labels.isna().any().any():
        raise ValueError("CellTypist predictions do not align with the input cell barcodes")
    confidence = pd.to_numeric(labels["conf_score"], errors="coerce").fillna(0.0)
    majority = labels["majority_voting"].astype(str)
    adata.obs["celltypist_type"] = labels["predicted_labels"].astype(str).to_numpy()
    adata.obs["celltypist_confidence"] = confidence.to_numpy()
    adata.obs["celltypist_cell_type"] = np.where(
        confidence.to_numpy() >= min_confidence,
        majority.to_numpy(),
        "Unknown",
    )


def load_manual(adata, tsv_path, cluster_key):
    """Load a complete, one-row-per-cluster reviewed mapping."""
    mapping = pd.read_csv(tsv_path, sep="\t", dtype=str)
    if not {"cluster", "cell_type"} <= set(mapping.columns):
        raise ValueError("manual TSV must have cluster and cell_type columns")
    mapping["cluster"] = mapping["cluster"].str.strip()
    mapping["cell_type"] = mapping["cell_type"].str.strip()
    if (
        mapping["cluster"].duplicated().any()
        or mapping[["cluster", "cell_type"]].isna().any().any()
        or mapping["cluster"].eq("").any()
        or mapping["cell_type"].eq("").any()
    ):
        raise ValueError("manual TSV clusters must be unique and cell_type cannot be blank")
    mapping = mapping.set_index("cluster")["cell_type"]
    clusters = set(adata.obs[cluster_key].astype(str).unique())
    missing = sorted(clusters - set(mapping.index.astype(str)))
    extra = sorted(set(mapping.index.astype(str)) - clusters)
    if missing:
        raise ValueError("manual TSV is missing clusters: " + ", ".join(missing))
    if extra:
        raise ValueError("manual TSV contains unknown clusters: " + ", ".join(extra))
    adata.obs["manual_cell_type"] = (
        adata.obs[cluster_key].astype(str).map(mapping.astype(str)).astype("category")
    )


def composition_table(adata, group_key, split_key):
    counts = pd.crosstab(adata.obs[split_key], adata.obs[group_key])
    proportions = counts.div(counts.sum(axis=1), axis=0)
    long_counts = counts.stack().rename("n_cells")
    long_proportions = proportions.stack().rename("proportion")
    return pd.concat([long_counts, long_proportions], axis=1).reset_index()


def composition_barplot(table, group_key, split_key, path):
    pivot = table.pivot(index=split_key, columns=group_key, values="proportion")
    fig, ax = plt.subplots(figsize=(max(8, len(pivot.index)), 5))
    pivot.plot.bar(stacked=True, ax=ax, colormap="tab20", width=0.8)
    ax.set_ylabel("Proportion")
    ax.set_xlabel(split_key)
    ax.set_title(f"Descriptive cell-type composition per {split_key}")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    adata = load_adata(args.input)
    if args.cluster_key not in adata.obs:
        raise ValueError(f"cluster key '{args.cluster_key}' is missing")

    plot_umap_grid(adata, [args.cluster_key], os.path.join(fig_dir, "01_umap_leiden.png"))
    marker_summary = pd.DataFrame()
    coverage = pd.DataFrame(
        columns=["cell_type", "markers_present", "markers_total", "eligible"]
    )
    marker_dictionary = None
    if args.markers_json:
        marker_dictionary = load_marker_dictionary(args.markers_json)
    elif args.marker_set:
        marker_dictionary = MARKER_SETS[args.marker_set]

    if marker_dictionary:
        score_keys, marker_summary, coverage = score_markers(
            adata,
            marker_dictionary,
            args.cluster_key,
            args.min_marker_coverage,
            args.min_score_margin,
            args.seed,
        )
        plot_umap_grid(
            adata,
            [key for _, key in score_keys[:6]],
            os.path.join(fig_dir, "02_umap_marker_scores.png"),
            title="Provisional marker-module scores",
        )
        plot_umap_grid(
            adata,
            ["marker_cell_type"],
            os.path.join(fig_dir, "03_umap_marker_cell_type.png"),
        )

    if args.celltypist:
        run_celltypist(adata, args.celltypist_model, args.celltypist_min_confidence)
        plot_umap_grid(
            adata,
            ["celltypist_cell_type", "celltypist_confidence"],
            os.path.join(fig_dir, "04_umap_celltypist.png"),
        )

    if args.manual_tsv:
        load_manual(adata, args.manual_tsv, args.cluster_key)
        plot_umap_grid(
            adata,
            ["manual_cell_type"],
            os.path.join(fig_dir, "05_umap_manual.png"),
        )

    if args.manual_tsv:
        source_column, method = "manual_cell_type", "manual"
    elif args.celltypist:
        source_column, method = "celltypist_cell_type", "celltypist"
    elif marker_dictionary:
        source_column, method = "marker_cell_type", "marker_score"
    else:
        source_column, method = None, "unassigned"
    adata.obs["cell_type"] = (
        adata.obs[source_column].astype(str) if source_column else pd.Series("Unknown", index=adata.obs_names)
    )
    adata.obs["annotation_method"] = method

    composition = composition_table(adata, "cell_type", "sample")
    composition.to_csv(os.path.join(args.out, "composition_by_sample.csv"), index=False)
    composition_barplot(
        composition,
        "cell_type",
        "sample",
        os.path.join(fig_dir, "06_composition_barplot.png"),
    )

    final_summary = (
        adata.obs.groupby([args.cluster_key, "cell_type"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    if not marker_summary.empty:
        final_summary[args.cluster_key] = final_summary[args.cluster_key].astype(str)
        marker_summary["cluster"] = marker_summary["cluster"].astype(str)
        final_summary = final_summary.merge(
            marker_summary,
            left_on=args.cluster_key,
            right_on="cluster",
            how="left",
        )
    final_summary["annotation_method"] = method
    final_summary.to_csv(os.path.join(args.out, "annotation_summary.csv"), index=False)
    coverage.to_csv(os.path.join(args.out, "marker_coverage.csv"), index=False)
    save_adata(adata, os.path.join(args.out, "annotated.h5ad"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--out", default="results/merged/08_annotate")
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--marker-set", choices=["", *MARKER_SETS], default="")
    parser.add_argument("--markers-json", default="")
    parser.add_argument("--min-marker-coverage", type=float, default=0.5)
    parser.add_argument("--min-score-margin", type=float, default=0.1)
    parser.add_argument("--celltypist", action="store_true")
    parser.add_argument("--celltypist-model", default="Immune_All_Low.pkl")
    parser.add_argument("--celltypist-min-confidence", type=float, default=0.5)
    parser.add_argument("--manual-tsv", default="")
    parser.add_argument("--seed", type=int, default=42)
    main(parser.parse_args())
