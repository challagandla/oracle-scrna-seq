"""
Core QC calculations: metric annotation, MAD filtering, stats.
"""
import json
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation


# ── Metric calculation ────────────────────────────────────────────────────────

def annotate_qc_vars(adata, mt_prefix="MT-", ribo_prefix="RPL,RPS",
                     hb_pattern=r"^HB[^(P)]"):
    """Flag mitochondrial, ribosomal, and hemoglobin genes in adata.var."""
    mt_tup   = tuple(mt_prefix.split(","))
    ribo_tup = tuple(ribo_prefix.split(","))
    adata.var["mt"]   = adata.var_names.str.startswith(mt_tup)
    adata.var["ribo"] = adata.var_names.str.startswith(ribo_tup)
    adata.var["hb"]   = adata.var_names.str.match(hb_pattern)
    return {
        "n_mt":   int(adata.var["mt"].sum()),
        "n_ribo": int(adata.var["ribo"].sum()),
        "n_hb":   int(adata.var["hb"].sum()),
    }


def calculate_qc_metrics(adata, **kwargs):
    """Wrapper around sc.pp.calculate_qc_metrics."""
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=[20, 50, 100, 200],  # pct_counts_in_top_N_genes
        log1p=False,
        inplace=True,
    )


# ── MAD filtering ─────────────────────────────────────────────────────────────

def mad_bounds(series, n_mads, log_transform=False):
    """Return (lower, upper) MAD bounds."""
    vals = np.log1p(series.values) if log_transform else series.values
    med  = np.median(vals)
    mad  = median_abs_deviation(vals)
    lower = med - n_mads * mad
    upper = med + n_mads * mad
    if log_transform:
        lower = np.expm1(lower)
        upper = np.expm1(upper)
    return lower, upper


def build_filter_mask(adata, qc_cfg):
    """
    Build a boolean pass_qc mask from config thresholds.
    Returns (pass_mask, thresholds_dict).
    """
    thresholds = {}

    # Total counts (log-transformed)
    lo_cnt, hi_cnt = mad_bounds(
        adata.obs["total_counts"], qc_cfg["mad_counts"], log_transform=True)
    out_counts = (
        (adata.obs["total_counts"] < max(lo_cnt, qc_cfg["min_counts"])) |
        (adata.obs["total_counts"] > hi_cnt)
    )
    thresholds["total_counts"] = {
        "lower": float(max(lo_cnt, qc_cfg["min_counts"])),
        "upper": float(hi_cnt),
    }

    # Gene counts (log-transformed)
    lo_g, hi_g = mad_bounds(
        adata.obs["n_genes_by_counts"], qc_cfg["mad_genes"], log_transform=True)
    out_genes = (
        (adata.obs["n_genes_by_counts"] < max(lo_g, qc_cfg["min_genes"])) |
        (adata.obs["n_genes_by_counts"] > hi_g)
    )
    thresholds["n_genes_by_counts"] = {
        "lower": float(max(lo_g, qc_cfg["min_genes"])),
        "upper": float(hi_g),
    }

    # MT % (raw, not log-transformed)
    lo_mt, hi_mt = mad_bounds(
        adata.obs["pct_counts_mt"], qc_cfg["mad_mt"], log_transform=False)
    hard_mt = float(qc_cfg["mt_hard"])
    effective_mt_upper = min(hi_mt, hard_mt)
    out_mt = adata.obs["pct_counts_mt"] > effective_mt_upper
    thresholds["pct_counts_mt"] = {
        "upper": float(effective_mt_upper),
        "hard":  hard_mt,
    }

    fail = out_counts | out_genes | out_mt
    return ~fail, thresholds


def summarise_qc(adata, label=""):
    stats = {
        "label":        label,
        "n_cells":      int(adata.n_obs),
        "n_genes":      int(adata.n_vars),
        "mean_counts":  float(adata.obs["total_counts"].mean()),
        "median_counts":float(adata.obs["total_counts"].median()),
        "mean_genes":   float(adata.obs["n_genes_by_counts"].mean()),
        "median_genes": float(adata.obs["n_genes_by_counts"].median()),
        "mean_pct_mt":  float(adata.obs["pct_counts_mt"].mean()),
    }
    print(f"\n── {label} ──")
    for k, v in stats.items():
        if k != "label":
            print(f"  {k:<20} {v:,.2f}" if isinstance(v, float) else f"  {k:<20} {v:,}")
    return stats
