#!/usr/bin/env python3
"""Cluster-marker ranking and replicate-aware pseudobulk condition DE."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from utils.io_utils import load_adata, save_adata
from utils.plot_utils import plot_dotplot, plot_heatmap, plot_stacked_violin

sc.settings.verbosity = 1
FILTERED_KEY = "rank_genes_groups_filtered"


def find_markers(adata, cluster_key, method, n_genes, logfc_min, pct_min):
    """Rank descriptive cluster markers and apply the configured log2FC filter."""
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method=method,
        use_raw=adata.raw is not None,
        pts=True,
        tie_correct=True,
        n_genes=n_genes,
    )
    sc.tl.filter_rank_genes_groups(
        adata,
        min_in_group_fraction=pct_min,
        min_fold_change=logfc_min,
        max_out_group_fraction=0.5,
        key_added=FILTERED_KEY,
    )
    adata.uns["marker_filter"] = {
        "min_in_group_fraction": float(pct_min),
        "min_log2_fold_change": float(logfc_min),
        "max_out_group_fraction": 0.5,
    }
    print(f"Marker detection complete for {adata.obs[cluster_key].nunique()} clusters")


def markers_to_df(adata, n_genes, key=FILTERED_KEY):
    """Convert a Scanpy ranking to a stable, H5AD-independent table."""
    ranking = adata.uns[key]
    groups = list(ranking["names"].dtype.names)
    full = adata.uns["rank_genes_groups"]
    fractions = ranking.get("pts", full.get("pts"))
    rows = []
    for group in groups:
        for index in range(min(n_genes, len(ranking["names"][group]))):
            gene = ranking["names"][group][index]
            if pd.isna(gene) or not str(gene).strip():
                continue
            fraction = np.nan
            if (
                isinstance(fractions, pd.DataFrame)
                and gene in fractions.index
                and group in fractions.columns
            ):
                fraction = fractions.loc[gene, group]
            rows.append(
                {
                    "cluster": group,
                    "gene": str(gene),
                    "score": ranking["scores"][group][index],
                    "log2FC": ranking["logfoldchanges"][group][index],
                    "padj": ranking["pvals_adj"][group][index],
                    "pct_grp": fraction,
                }
            )
    return pd.DataFrame(rows, columns=["cluster", "gene", "score", "log2FC", "padj", "pct_grp"])


def volcano_plot(table, cluster, path):
    """Plot descriptive one-vs-rest cluster-marker statistics."""
    fig, ax = plt.subplots(figsize=(5, 4))
    x_values = table["log2FC"].fillna(0)
    y_values = (-np.log10(table["padj"].replace(0, 1e-300))).fillna(0)
    colours = np.where(
        (x_values > 0.25) & (table["padj"] < 0.05),
        "red",
        np.where((x_values < -0.25) & (table["padj"] < 0.05), "blue", "grey"),
    )
    ax.scatter(x_values, y_values, c=colours, s=5, alpha=0.6, rasterized=True)
    ax.axhline(-np.log10(0.05), ls="--", color="black", lw=0.8)
    ax.axvline(0.25, ls="--", color="black", lw=0.8)
    ax.axvline(-0.25, ls="--", color="black", lw=0.8)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(f"Cluster {cluster} vs all other cells (descriptive)")
    for _, row in table.nsmallest(10, "padj").iterrows():
        fold_change = row["log2FC"] if np.isfinite(row["log2FC"]) else 0
        p_value = row["padj"] if np.isfinite(row["padj"]) else 1
        ax.annotate(
            row["gene"],
            (fold_change, -np.log10(max(p_value, 1e-300))),
            fontsize=6,
            ha="center",
        )
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _constant_sample_value(obs, sample, column):
    values = obs.loc[obs["sample"].astype(str) == str(sample), column].dropna().astype(str).unique()
    if len(values) != 1:
        raise ValueError(f"sample '{sample}' must have exactly one value for '{column}'")
    return values[0]


def pseudobulk_de(
    adata,
    cluster_key,
    condition_key,
    group1,
    group2,
    de_dir,
    *,
    covariates,
    min_cells_per_sample,
    min_replicates_per_group,
    min_total_count,
):
    """Run per-cluster sample-level DE and return an explicit status table."""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as error:
        raise RuntimeError("PyDESeq2 is required when condition DE is configured") from error

    required = {"sample", condition_key, *covariates}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise ValueError(f"pseudobulk metadata columns are missing: {', '.join(missing)}")
    if "counts" not in adata.layers:
        raise ValueError("pseudobulk DE requires the preserved raw-count layer")

    os.makedirs(de_dir, exist_ok=True)
    statuses = []
    design_terms = [*covariates, condition_key]
    design = "~" + " + ".join(design_terms)

    for cluster in sorted(adata.obs[cluster_key].astype(str).unique()):
        cluster_data = adata[
            (adata.obs[cluster_key].astype(str) == cluster)
            & adata.obs[condition_key].astype(str).isin([group1, group2])
        ]
        count_rows = []
        metadata_rows = []
        excluded_samples = []
        eligible_samples = adata.obs.loc[
            adata.obs[condition_key].astype(str).isin([group1, group2]), "sample"
        ].astype(str).unique()
        for sample in sorted(eligible_samples):
            sample_data = cluster_data[cluster_data.obs["sample"].astype(str) == sample]
            if sample_data.n_obs < min_cells_per_sample:
                excluded_samples.append(sample)
                continue
            matrix = sample_data.layers["counts"]
            aggregate = (
                np.asarray(matrix.sum(axis=0)).ravel()
                if sparse.issparse(matrix)
                else np.asarray(matrix).sum(axis=0)
            )
            if not np.isfinite(aggregate).all() or (aggregate < 0).any() or not np.allclose(
                aggregate, np.rint(aggregate), atol=1e-6, rtol=0
            ):
                raise ValueError(f"cluster {cluster}: pseudobulk counts are not non-negative integers")
            metadata = {"sample": sample}
            for column in design_terms:
                metadata[column] = _constant_sample_value(adata.obs, sample, column)
            count_rows.append(np.rint(aggregate).astype(np.int64))
            metadata_rows.append(metadata)

        metadata = pd.DataFrame(metadata_rows)
        replicate_counts = (
            metadata[condition_key].value_counts() if not metadata.empty else pd.Series(dtype=int)
        )
        reason = ""
        for group in (group1, group2):
            if int(replicate_counts.get(group, 0)) < min_replicates_per_group:
                reason = (
                    f"fewer than {min_replicates_per_group} usable biological replicates "
                    f"for group '{group}'"
                )
                break
        if reason:
            statuses.append(
                {
                    "cluster": cluster,
                    "status": "skipped",
                    "reason": reason,
                    "n_pseudobulks": len(metadata_rows),
                    "excluded_samples": ",".join(excluded_samples),
                    "output": "",
                }
            )
            continue

        counts = pd.DataFrame(
            np.vstack(count_rows),
            index=metadata["sample"],
            columns=cluster_data.var_names,
        )
        keep_genes = counts.sum(axis=0) >= min_total_count
        counts = counts.loc[:, keep_genes]
        metadata = metadata.set_index("sample")
        if counts.shape[1] < 2:
            statuses.append(
                {
                    "cluster": cluster,
                    "status": "skipped",
                    "reason": "fewer than two genes passed the pseudobulk count filter",
                    "n_pseudobulks": len(metadata),
                    "excluded_samples": ",".join(excluded_samples),
                    "output": "",
                }
            )
            continue

        try:
            try:
                dataset = DeseqDataSet(counts=counts, metadata=metadata, design=design, quiet=True)
            except TypeError:
                dataset = DeseqDataSet(
                    counts=counts,
                    metadata=metadata,
                    design_factors=design_terms,
                    quiet=True,
                )
            dataset.deseq2()
            statistics = DeseqStats(
                dataset,
                contrast=[condition_key, group2, group1],
                quiet=True,
            )
            statistics.summary()
            results = statistics.results_df.sort_values("padj")
            safe_cluster = re.sub(r"[^A-Za-z0-9_.-]+", "_", cluster)
            output = os.path.join(de_dir, f"pseudobulk_cluster_{safe_cluster}.csv")
            results.to_csv(output)
            statuses.append(
                {
                    "cluster": cluster,
                    "status": "complete",
                    "reason": f"log2FC is {group2} versus reference {group1}",
                    "n_pseudobulks": len(metadata),
                    "excluded_samples": ",".join(excluded_samples),
                    "output": output,
                }
            )
        except Exception as error:  # preserve a visible per-cluster audit trail
            statuses.append(
                {
                    "cluster": cluster,
                    "status": "failed",
                    "reason": str(error),
                    "n_pseudobulks": len(metadata),
                    "excluded_samples": ",".join(excluded_samples),
                    "output": "",
                }
            )

    return pd.DataFrame(statuses)


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    de_dir = os.path.join(args.out, "DE")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(de_dir, exist_ok=True)
    for stale_result in Path(de_dir).glob("pseudobulk_cluster_*.csv"):
        stale_result.unlink()
    adata = load_adata(args.input)
    if args.cluster_key not in adata.obs:
        raise ValueError(f"cluster key '{args.cluster_key}' is missing")

    find_markers(
        adata,
        args.cluster_key,
        method=args.method,
        n_genes=args.n_genes,
        logfc_min=args.logfc_min,
        pct_min=args.pct_min,
    )
    marker_table = markers_to_df(adata, args.n_genes)
    marker_table.to_csv(os.path.join(args.out, "cluster_markers.csv"), index=False)

    for plotting_function, filename, n_genes in (
        (plot_dotplot, "01_markers_dotplot.png", 5),
        (plot_heatmap, "02_markers_heatmap.png", 5),
        (plot_stacked_violin, "03_markers_stacked_violin.png", 3),
    ):
        try:
            plotting_function(
                adata,
                args.cluster_key,
                n_genes=n_genes,
                path=os.path.join(fig_dir, filename),
            )
        except Exception as error:
            print(f"{plotting_function.__name__} failed: {error}")

    for cluster in adata.obs[args.cluster_key].value_counts().head(5).index:
        subset = marker_table[marker_table["cluster"].astype(str) == str(cluster)]
        if len(subset) > 1:
            volcano_plot(subset, cluster, os.path.join(fig_dir, f"04_volcano_cluster_{cluster}.png"))

    status_path = os.path.join(args.out, "pseudobulk_status.csv")
    if args.condition_key:
        try:
            status = pseudobulk_de(
                adata,
                args.cluster_key,
                args.condition_key,
                args.group1,
                args.group2,
                de_dir,
                covariates=args.covariates,
                min_cells_per_sample=args.min_cells_per_sample,
                min_replicates_per_group=args.min_replicates_per_group,
                min_total_count=args.min_total_count,
            )
        except Exception as error:
            fatal_status = pd.DataFrame(
                [{"cluster": "", "status": "failed", "reason": str(error)}]
            )
            fatal_status.to_csv(status_path, index=False)
            print("Pseudobulk setup failed:")
            print(fatal_status.to_csv(index=False))
            raise
    else:
        status = pd.DataFrame(
            [{"cluster": "", "status": "not_requested", "reason": "condition_key is empty"}]
        )
    status.to_csv(status_path, index=False)
    if args.condition_key and not status["status"].eq("complete").any():
        print("Pseudobulk status (also shown here because failed-job outputs may be removed):")
        print(status.to_csv(index=False))
        details = "; ".join(
            f"cluster {row.cluster}: {row.reason}"
            for row in status.head(5).itertuples(index=False)
        )
        raise RuntimeError(f"No pseudobulk comparison completed. {details}")

    # The filtered Scanpy structure may mix NaN and string placeholders. Keep
    # the canonical ranking plus an explicit filter contract in the H5AD, and
    # store the filtered values in the CSV written above.
    adata.uns.pop(FILTERED_KEY, None)
    save_adata(adata, os.path.join(args.out, "with_markers.h5ad"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--out", default="results/merged/07_markers")
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--method", default="wilcoxon", choices=["wilcoxon", "t-test"])
    parser.add_argument("--n-genes", type=int, default=50)
    parser.add_argument("--logfc-min", type=float, default=0.25)
    parser.add_argument("--pct-min", type=float, default=0.1)
    parser.add_argument("--condition-key", default="")
    parser.add_argument("--group1", default="")
    parser.add_argument("--group2", default="")
    parser.add_argument("--covariates", nargs="*", default=[])
    parser.add_argument("--min-cells-per-sample", type=int, default=10)
    parser.add_argument("--min-replicates-per-group", type=int, default=2)
    parser.add_argument("--min-total-count", type=int, default=10)
    main(parser.parse_args())
