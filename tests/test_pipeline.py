"""Lightweight unit tests for pipeline utility functions."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from utils.qc_utils import (
    annotate_qc_vars,
    build_filter_mask,
    calculate_qc_metrics,
    mad_bounds,
    summarise_qc,
)


# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def tiny_adata():
    """200 cells × 100 genes with synthetic counts."""
    np.random.seed(42)
    n_cells, n_genes = 200, 100
    counts = np.random.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(float)
    gene_names = [f"Gene{i:04d}" for i in range(n_genes)]
    # Add a handful of MT genes
    gene_names[0] = "MT-ND1"
    gene_names[1] = "MT-CO1"
    gene_names[2] = "MT-ATP6"
    adata = ad.AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame(index=[f"Cell{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=gene_names),
    )
    return adata


# ── QC utils ─────────────────────────────────────────────────────────────────

def test_annotate_qc_vars(tiny_adata):
    result = annotate_qc_vars(tiny_adata, mt_prefix="MT-")
    assert "mt" in tiny_adata.var.columns
    assert "ribo" in tiny_adata.var.columns
    assert "hb" in tiny_adata.var.columns
    assert result["n_mt"] == 3


def test_calculate_qc_metrics(tiny_adata):
    annotate_qc_vars(tiny_adata, mt_prefix="MT-")
    calculate_qc_metrics(tiny_adata)
    assert "total_counts" in tiny_adata.obs.columns
    assert "pct_counts_mt" in tiny_adata.obs.columns
    assert tiny_adata.obs["total_counts"].min() > 0


@pytest.mark.parametrize(
    ("n_vars", "expected"),
    [(19, None), (21, [20]), (75, [20, 50]), (101, [20, 50, 100])],
)
def test_calculate_qc_metrics_bounds_percent_top(monkeypatch, n_vars, expected):
    captured = {}

    class TinyShape:
        pass

    adata = TinyShape()
    adata.n_vars = n_vars

    def fake_calculate(*args, **kwargs):
        captured["percent_top"] = kwargs["percent_top"]

    monkeypatch.setattr("utils.qc_utils.sc.pp.calculate_qc_metrics", fake_calculate)
    calculate_qc_metrics(adata)
    assert captured["percent_top"] == expected


def test_mad_bounds_log(tiny_adata):
    import scanpy as sc
    annotate_qc_vars(tiny_adata, mt_prefix="MT-")
    sc.pp.calculate_qc_metrics(
        tiny_adata, qc_vars=["mt", "ribo", "hb"],
        percent_top=None, log1p=False, inplace=True)
    lo, hi = mad_bounds(tiny_adata.obs["total_counts"], n_mads=5, log_transform=True)
    assert lo < hi
    assert lo > 0


def test_build_filter_mask_keeps_most_cells(tiny_adata):
    import scanpy as sc
    annotate_qc_vars(tiny_adata, mt_prefix="MT-")
    sc.pp.calculate_qc_metrics(
        tiny_adata, qc_vars=["mt", "ribo", "hb"],
        percent_top=None, log1p=False, inplace=True)
    qc_cfg = {
        "mad_counts": 5, "mad_genes": 5, "mad_mt": 3,
        "mt_hard": 80, "min_counts": 0, "min_genes": 0,
    }
    mask, thresholds = build_filter_mask(tiny_adata, qc_cfg)
    # With permissive thresholds, >90% cells should pass
    assert mask.mean() > 0.9
    assert "total_counts" in thresholds
    assert "pct_counts_mt" in thresholds


def test_build_filter_mask_stringent(tiny_adata):
    import scanpy as sc
    annotate_qc_vars(tiny_adata, mt_prefix="MT-")
    sc.pp.calculate_qc_metrics(
        tiny_adata, qc_vars=["mt", "ribo", "hb"],
        percent_top=None, log1p=False, inplace=True)
    qc_cfg = {
        "mad_counts": 0.1, "mad_genes": 0.1, "mad_mt": 0.1,
        "mt_hard": 0.1, "min_counts": 9999, "min_genes": 9999,
    }
    mask, _ = build_filter_mask(tiny_adata, qc_cfg)
    # Very stringent — should remove most cells
    assert mask.mean() < 0.5


def test_summarise_qc(tiny_adata, capsys):
    import scanpy as sc
    annotate_qc_vars(tiny_adata, mt_prefix="MT-")
    sc.pp.calculate_qc_metrics(
        tiny_adata, qc_vars=["mt", "ribo", "hb"],
        percent_top=None, log1p=False, inplace=True)
    stats = summarise_qc(tiny_adata, label="Test")
    assert stats["n_cells"] == 200
    assert stats["n_genes"] == 100
    assert "mean_counts" in stats


# ── I/O ──────────────────────────────────────────────────────────────────────

def test_load_save_roundtrip(tmp_path, tiny_adata):
    from utils.io_utils import save_adata, load_adata
    path = str(tmp_path / "test.h5ad")
    save_adata(tiny_adata, path)
    loaded = load_adata(path)
    assert loaded.n_obs == tiny_adata.n_obs
    assert loaded.n_vars == tiny_adata.n_vars


def test_save_adata_keeps_existing_file_when_write_fails(tmp_path, tiny_adata, monkeypatch):
    from utils.io_utils import save_adata

    path = tmp_path / "stable.h5ad"
    path.write_bytes(b"previous-complete-result")

    def fail_after_partial_write(self, destination):
        del self
        Path(destination).write_bytes(b"partial")
        raise RuntimeError("simulated interrupted write")

    monkeypatch.setattr(ad.AnnData, "write_h5ad", fail_after_partial_write)
    with pytest.raises(RuntimeError, match="simulated interrupted write"):
        save_adata(tiny_adata, str(path))

    assert path.read_bytes() == b"previous-complete-result"
    assert not list(tmp_path.glob(".*.tmp.h5ad"))
