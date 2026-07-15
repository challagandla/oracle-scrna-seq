"""
Reusable plotting helpers for the scRNA-seq pipeline.
All functions save to disk and return the figure path.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

PALETTE = "tab20"


def _savefig(fig, path, dpi=150):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


# ── QC ───────────────────────────────────────────────────────────────────────

def plot_qc_violin(adata, metrics, title, path):
    """Violin plots for QC metrics."""
    n = len(metrics)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4))
    if n == 1:
        axes = [axes]
    for ax, m in zip(axes, metrics, strict=True):
        data = adata.obs[m].dropna()
        ax.violinplot(data, showmedians=True)
        ax.set_title(m.replace("_", " "), fontsize=9)
        ax.set_xticks([])
    fig.suptitle(title, fontsize=11, y=1.02)
    plt.tight_layout()
    return _savefig(fig, path)


def plot_qc_scatter(adata, x, y, color, title, path):
    """Scatter plot of two QC metrics coloured by a third."""
    fig, ax = plt.subplots(figsize=(5, 4))
    vals = adata.obs[color]
    sc_plot = ax.scatter(
        adata.obs[x], adata.obs[y],
        c=vals, cmap="RdYlBu_r", s=1, alpha=0.4, rasterized=True
    )
    plt.colorbar(sc_plot, ax=ax, label=color.replace("_", " "))
    ax.set_xlabel(x.replace("_", " "))
    ax.set_ylabel(y.replace("_", " "))
    ax.set_title(title)
    return _savefig(fig, path)


def plot_qc_histograms(adata, metrics, thresholds, path):
    """
    Histogram per metric with threshold lines overlaid.
    thresholds: dict  metric -> {"lower": val, "upper": val, "hard": val}
    """
    n = len(metrics)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 3.5))
    if n == 1:
        axes = [axes]
    for ax, m in zip(axes, metrics, strict=True):
        vals = adata.obs[m].dropna()
        ax.hist(vals, bins=60, color="steelblue", edgecolor="none", alpha=0.7)
        t = thresholds.get(m, {})
        if "lower" in t:
            ax.axvline(t["lower"], color="red",    ls="--", lw=1.2, label=f'lower={t["lower"]:.1f}')
        if "upper" in t:
            ax.axvline(t["upper"], color="orange", ls="--", lw=1.2, label=f'upper={t["upper"]:.1f}')
        if "hard" in t:
            ax.axvline(t["hard"],  color="black",  ls="-",  lw=1.2, label=f'hard={t["hard"]:.1f}')
        ax.set_xlabel(m.replace("_", " "), fontsize=9)
        ax.set_ylabel("Cells")
        ax.legend(fontsize=7)
    plt.tight_layout()
    return _savefig(fig, path)


def plot_knee(adata, path):
    """Knee/barcode-rank plot (like Cell Ranger's web summary)."""
    total_counts = np.sort(adata.obs["total_counts"].values)[::-1]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.loglog(np.arange(1, len(total_counts) + 1), total_counts, lw=1.5)
    ax.set_xlabel("Barcode rank")
    ax.set_ylabel("Total UMI counts")
    ax.set_title("Knee plot")
    ax.grid(True, which="both", alpha=0.3)
    return _savefig(fig, path)


# ── HVG / PCA ────────────────────────────────────────────────────────────────

def plot_elbow(adata, n_pcs_use, path):
    """PCA elbow plot with cumulative variance and recommended cutoff."""
    vr  = adata.uns["pca"]["variance_ratio"]
    cum = np.cumsum(vr) * 100
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    axes[0].plot(range(1, len(vr) + 1), vr, "o-", ms=3, lw=1)
    axes[0].axvline(n_pcs_use, color="red", ls="--", label=f"n_pcs={n_pcs_use}")
    axes[0].set_xlabel("PC")
    axes[0].set_ylabel("Variance ratio")
    axes[0].set_title("Elbow plot")
    axes[0].legend()

    axes[1].plot(range(1, len(cum) + 1), cum, "o-", ms=3, lw=1)
    axes[1].axhline(80, color="orange", ls="--", alpha=0.7, label="80%")
    axes[1].axhline(90, color="red",    ls="--", alpha=0.7, label="90%")
    axes[1].axvline(n_pcs_use, color="blue", ls="--", label=f"n_pcs={n_pcs_use}")
    axes[1].set_xlabel("PC")
    axes[1].set_ylabel("Cumulative variance (%)")
    axes[1].set_title("Cumulative variance")
    axes[1].legend(fontsize=7)

    plt.tight_layout()
    return _savefig(fig, path)


# ── UMAP ─────────────────────────────────────────────────────────────────────

def plot_umap_grid(adata, color_keys, path, ncols=3, title=""):
    """Multi-panel UMAP grid coloured by obs columns or gene names."""
    valid = [k for k in color_keys
             if k in adata.obs.columns or k in adata.var_names]
    if not valid:
        return None
    nrows = (len(valid) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(5 * ncols, 4 * nrows), squeeze=False)
    flat = [ax for row in axes for ax in row]
    for i, k in enumerate(valid):
        sc.pl.umap(adata, color=k, ax=flat[i], show=False,
                   frameon=False, title=k, s=5)
    for j in range(i + 1, len(flat)):
        flat[j].set_visible(False)
    if title:
        fig.suptitle(title, y=1.01)
    plt.tight_layout()
    return _savefig(fig, path)


# ── Clustering ───────────────────────────────────────────────────────────────

def plot_silhouette(df, configured_res, path):
    valid = df["silhouette"].dropna()
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(valid.index, valid.values, "o-", lw=1.5)
    ax.axvline(
        configured_res,
        color="red",
        ls="--",
        label=f"configured={configured_res}",
    )
    ax.set_xlabel("Resolution")
    ax.set_ylabel("Mean silhouette score")
    ax.set_title("Leiden resolution diagnostics")
    ax.legend()
    plt.tight_layout()
    return _savefig(fig, path)


# ── Markers ──────────────────────────────────────────────────────────────────

def plot_dotplot(adata, cluster_key, n_genes, path):
    sc.pl.rank_genes_groups_dotplot(
        adata, groupby=cluster_key, n_genes=n_genes,
        show=False, return_fig=False,
    )
    fig = plt.gcf()
    return _savefig(fig, path)


def plot_heatmap(adata, cluster_key, n_genes, path):
    sc.pl.rank_genes_groups_heatmap(
        adata, groupby=cluster_key, n_genes=n_genes,
        show=False, save=False, swap_axes=True,
    )
    fig = plt.gcf()
    return _savefig(fig, path)


def plot_stacked_violin(adata, cluster_key, n_genes, path):
    sc.pl.rank_genes_groups_stacked_violin(
        adata, groupby=cluster_key, n_genes=n_genes,
        show=False, return_fig=False,
    )
    fig = plt.gcf()
    return _savefig(fig, path)
