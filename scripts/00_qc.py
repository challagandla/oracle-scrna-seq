#!/usr/bin/env python3
"""
Step 00 — Quality Control
=========================
Per-sample QC following scverse best practices.

Outputs
-------
results/<sample>/00_qc/
  qc_filtered.h5ad        clean AnnData
  qc_stats.json           before/after statistics
  figures/
    01_knee_plot.png       barcode-rank (knee) plot
    02_qc_violin_before.png
    03_qc_histograms.png   per-metric histograms with MAD thresholds
    04_counts_vs_genes.png scatter coloured by MT%
    05_counts_vs_mt.png
    06_qc_violin_after.png
"""

import argparse, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import scanpy as sc
import numpy as np

from utils.io_utils   import load_adata, save_adata
from utils.qc_utils   import (annotate_qc_vars, calculate_qc_metrics,
                               build_filter_mask, summarise_qc)
from utils.plot_utils import (plot_knee, plot_qc_violin, plot_qc_histograms,
                               plot_qc_scatter)

sc.settings.verbosity = 1


def main(args):
    fig_dir = os.path.join(args.out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    # ── Load ─────────────────────────────────────────────────
    adata = load_adata(args.input)

    # ── Knee plot ────────────────────────────────────────────
    # Needs total_counts — do a quick calculation first
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    plot_knee(adata, os.path.join(fig_dir, "01_knee_plot.png"))

    # ── Full QC metrics ──────────────────────────────────────
    gene_counts = annotate_qc_vars(
        adata,
        mt_prefix=args.mt_prefix,
        ribo_prefix=args.ribo_prefix,
        hb_pattern=args.hb_pattern,
    )
    print(f"  MT genes: {gene_counts['n_mt']}  |  "
          f"Ribo genes: {gene_counts['n_ribo']}  |  "
          f"Hb genes: {gene_counts['n_hb']}")
    calculate_qc_metrics(adata)

    qc_metrics = ["total_counts", "n_genes_by_counts",
                  "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]

    # ── Pre-filter plots ─────────────────────────────────────
    plot_qc_violin(adata, qc_metrics, "Before QC filtering",
                   os.path.join(fig_dir, "02_qc_violin_before.png"))

    # ── Build filter mask + threshold dict ───────────────────
    qc_cfg = {
        "mad_counts": args.mad_counts,
        "mad_genes":  args.mad_genes,
        "mad_mt":     args.mad_mt,
        "mt_hard":    args.mt_hard,
        "min_counts": args.min_counts,
        "min_genes":  args.min_genes,
    }
    pass_mask, thresholds = build_filter_mask(adata, qc_cfg)

    # Histogram with thresholds
    plot_qc_histograms(
        adata,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
        thresholds,
        os.path.join(fig_dir, "03_qc_histograms.png"),
    )

    # Scatter plots
    plot_qc_scatter(adata, "total_counts", "n_genes_by_counts",
                    "pct_counts_mt", "Counts vs Genes (colour=MT%)",
                    os.path.join(fig_dir, "04_counts_vs_genes.png"))
    plot_qc_scatter(adata, "total_counts", "pct_counts_mt",
                    "n_genes_by_counts", "Counts vs MT% (colour=genes)",
                    os.path.join(fig_dir, "05_counts_vs_mt.png"))

    # ── Before stats ─────────────────────────────────────────
    stats_before = summarise_qc(adata, label="Before QC")

    # ── Apply filters ────────────────────────────────────────
    n_fail = (~pass_mask).sum()
    print(f"\nCells failing QC: {n_fail} / {adata.n_obs} "
          f"({n_fail / adata.n_obs * 100:.1f}%)")

    adata.obs["pass_qc"] = pass_mask
    adata_f = adata[pass_mask].copy()

    # Gene filter
    sc.pp.filter_genes(adata_f, min_cells=args.min_cells)
    print(f"Genes after filter: {adata_f.n_vars:,}")

    stats_after = summarise_qc(adata_f, label="After QC")

    # ── Post-filter plots ────────────────────────────────────
    plot_qc_violin(adata_f, ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
                   "After QC filtering",
                   os.path.join(fig_dir, "06_qc_violin_after.png"))

    # ── Save stats JSON ──────────────────────────────────────
    stats = {
        "before": stats_before,
        "after":  stats_after,
        "thresholds": thresholds,
        "gene_counts": gene_counts,
    }
    with open(os.path.join(args.out, "qc_stats.json"), "w") as f:
        json.dump(stats, f, indent=2)

    save_adata(adata_f, os.path.join(args.out, "qc_filtered.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Per-sample QC")
    p.add_argument("input")
    p.add_argument("--out",          default="results/qc")
    p.add_argument("--mt-prefix",    default="MT-")
    p.add_argument("--ribo-prefix",  default="RPL,RPS")
    p.add_argument("--hb-pattern",   default=r"^HB[^(P)]")
    p.add_argument("--mad-counts",   type=float, default=5)
    p.add_argument("--mad-genes",    type=float, default=5)
    p.add_argument("--mad-mt",       type=float, default=3)
    p.add_argument("--mt-hard",      type=float, default=8)
    p.add_argument("--min-cells",    type=int,   default=20)
    p.add_argument("--min-counts",   type=int,   default=500)
    p.add_argument("--min-genes",    type=int,   default=200)
    main(p.parse_args())
