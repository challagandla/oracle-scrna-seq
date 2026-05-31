#!/usr/bin/env python3
"""
Step 07 — Marker Genes & Differential Expression
==================================================
A) Cluster markers via Wilcoxon one-vs-rest (recommended by scverse)
B) Condition DE via pseudobulk + pydeseq2 (requires biological replicates)

Outputs
-------
results/merged/07_markers/
  cluster_markers.csv         all markers per cluster
  with_markers.h5ad
  figures/
    01_markers_dotplot.png
    02_markers_heatmap.png
    03_markers_stacked_violin.png
    04_volcano_cluster_<N>.png  (top clusters)
  DE/
    pseudobulk_cluster_<N>.csv  (if condition DE requested)
"""

import argparse, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.io_utils   import load_adata, save_adata
from utils.plot_utils import plot_dotplot, plot_heatmap, plot_stacked_violin

sc.settings.verbosity = 1


def find_markers(adata, cluster_key, method, n_genes, logfc_min, pct_min):
    sc.tl.rank_genes_groups(
        adata, groupby=cluster_key, method=method,
        use_raw=(adata.raw is not None),
        pts=True, tie_correct=True,
        n_genes=n_genes,
    )
    # Apply logFC and pct filter
    sc.tl.filter_rank_genes_groups(
        adata,
        min_in_group_fraction=pct_min,
        min_fold_change=2 ** logfc_min,
        max_out_group_fraction=0.5,
    )
    print(f"Marker detection complete for {adata.obs[cluster_key].nunique()} clusters")


def markers_to_df(adata, n_genes):
    groups = list(adata.uns["rank_genes_groups"]["names"].dtype.names)
    rows   = []
    rgg    = adata.uns["rank_genes_groups"]
    for g in groups:
        for i in range(min(n_genes, len(rgg["names"][g]))):
            rows.append({
                "cluster":  g,
                "gene":     rgg["names"][g][i],
                "score":    rgg["scores"][g][i],
                "log2FC":   rgg["logfoldchanges"][g][i],
                "padj":     rgg["pvals_adj"][g][i],
                "pct_grp":  rgg.get("pts", {}).get(g, [np.nan]*n_genes)[i],
            })
    return pd.DataFrame(rows)


def volcano_plot(df_cluster, cluster, path):
    """Volcano plot for a single cluster's markers."""
    fig, ax = plt.subplots(figsize=(5, 4))
    x = df_cluster["log2FC"].fillna(0)
    y = (-np.log10(df_cluster["padj"].replace(0, 1e-300))).fillna(0)
    colors = np.where((x > 0.25) & (df_cluster["padj"] < 0.05), "red",
             np.where((x < -0.25) & (df_cluster["padj"] < 0.05), "blue", "grey"))
    ax.scatter(x, y, c=colors, s=5, alpha=0.6, rasterized=True)
    ax.axhline(-np.log10(0.05), ls="--", color="black", lw=0.8)
    ax.axvline(0.25,  ls="--", color="black", lw=0.8)
    ax.axvline(-0.25, ls="--", color="black", lw=0.8)
    ax.set_xlabel("log2 Fold Change"); ax.set_ylabel("-log10 adj. p-value")
    ax.set_title(f"Cluster {cluster} — volcano")
    # Label top 10
    top = df_cluster.nsmallest(10, "padj")
    for _, row in top.iterrows():
        lfc = row["log2FC"] if not np.isnan(row["log2FC"]) else 0
        pv  = row["padj"]   if not np.isnan(row["padj"])   else 1
        ax.annotate(row["gene"], (lfc, -np.log10(max(pv, 1e-300))),
                    fontsize=6, ha="center")
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def pseudobulk_de(adata, cluster_key, condition_key, group1, group2, de_dir):
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds  import DeseqStats
    except ImportError:
        print("pydeseq2 not installed — skipping pseudobulk DE. pip install pydeseq2")
        return

    if "sample" not in adata.obs.columns:
        print("No 'sample' column — pseudobulk requires biological replicates. Skipping.")
        return

    import scipy.sparse as sp
    os.makedirs(de_dir, exist_ok=True)

    for cluster in sorted(adata.obs[cluster_key].unique()):
        sub = adata[adata.obs[cluster_key] == cluster]
        sub = sub[sub.obs[condition_key].isin([group1, group2])]
        if sub.n_obs < 20:
            continue

        samples = sub.obs["sample"].unique()
        counts, meta_rows = [], []
        for s in samples:
            s_adata = sub[sub.obs["sample"] == s]
            raw = s_adata.layers.get("counts", s_adata.X)
            agg = np.asarray(raw.sum(axis=0)).flatten() if sp.issparse(raw) \
                  else np.asarray(raw).sum(axis=0)
            counts.append(agg)
            meta_rows.append({
                "sample": s,
                condition_key: s_adata.obs[condition_key].iloc[0]
            })

        counts_df = pd.DataFrame(
            np.array(counts),
            index=[r["sample"] for r in meta_rows],
            columns=sub.var_names,
        )
        meta_df = pd.DataFrame(meta_rows).set_index("sample")
        meta_df[condition_key] = meta_df[condition_key].astype("category")

        try:
            dds  = DeseqDataSet(counts=counts_df, metadata=meta_df,
                                design_factors=condition_key, quiet=True)
            dds.deseq2()
            stat = DeseqStats(dds, contrast=[condition_key, group1, group2], quiet=True)
            stat.summary()
            res  = stat.results_df.sort_values("padj")
            out  = os.path.join(de_dir, f"pseudobulk_cluster_{cluster}.csv")
            res.to_csv(out)
            sig  = (res["padj"] < 0.05).sum()
            print(f"  Cluster {cluster}: {sig} DE genes (padj<0.05) → {out}")
        except Exception as e:
            print(f"  Cluster {cluster}: DE failed — {e}")


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_adata(args.input)

    cluster_key = args.cluster_key
    if cluster_key not in adata.obs.columns:
        sys.exit(f"'{cluster_key}' not in adata.obs. Available: {list(adata.obs.columns)}")

    # ── Marker detection ─────────────────────────────────────
    find_markers(adata, cluster_key, method=args.method,
                 n_genes=args.n_genes, logfc_min=args.logfc_min,
                 pct_min=args.pct_min)

    df = markers_to_df(adata, n_genes=args.n_genes)
    df.to_csv(os.path.join(args.out, "cluster_markers.csv"), index=False)

    # ── Visualizations ────────────────────────────────────────
    try:
        plot_dotplot(adata, cluster_key, n_genes=5,
                     path=os.path.join(fig_dir, "01_markers_dotplot.png"))
    except Exception as e:
        print(f"Dotplot failed: {e}")

    try:
        plot_heatmap(adata, cluster_key, n_genes=5,
                     path=os.path.join(fig_dir, "02_markers_heatmap.png"))
    except Exception as e:
        print(f"Heatmap failed: {e}")

    try:
        plot_stacked_violin(adata, cluster_key, n_genes=3,
                            path=os.path.join(fig_dir, "03_markers_stacked_violin.png"))
    except Exception as e:
        print(f"Stacked violin failed: {e}")

    # Volcano per cluster (top 5 clusters by size)
    top_clusters = (adata.obs[cluster_key].value_counts()
                    .head(5).index.tolist())
    for cl in top_clusters:
        sub = df[df["cluster"] == cl]
        if len(sub) > 1:
            volcano_plot(sub, cl,
                         os.path.join(fig_dir, f"04_volcano_cluster_{cl}.png"))

    # ── Condition DE ──────────────────────────────────────────
    if args.condition_key and args.group1 and args.group2:
        print(f"\nRunning pseudobulk DE: {args.group1} vs {args.group2}")
        pseudobulk_de(adata, cluster_key, args.condition_key,
                      args.group1, args.group2,
                      os.path.join(args.out, "DE"))

    save_adata(adata, os.path.join(args.out, "with_markers.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--out",           default="results/merged/07_markers")
    p.add_argument("--cluster-key",   default="leiden")
    p.add_argument("--method",        default="wilcoxon",
                   choices=["wilcoxon", "t-test", "logreg"])
    p.add_argument("--n-genes",       type=int,   default=50)
    p.add_argument("--logfc-min",     type=float, default=0.25)
    p.add_argument("--pct-min",       type=float, default=0.1)
    p.add_argument("--condition-key", default="")
    p.add_argument("--group1",        default="")
    p.add_argument("--group2",        default="")
    main(p.parse_args())
