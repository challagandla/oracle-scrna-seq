"""Focused regressions for marker, annotation, and trajectory contracts."""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


def clustered_adata():
    return ad.AnnData(
        X=sparse.csr_matrix(np.ones((4, 4), dtype=np.float32)),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(["0", "0", "1", "1"])},
            index=[f"cell_{index}" for index in range(4)],
        ),
        var=pd.DataFrame(index=["A1", "A2", "B1", "B2"]),
    )


def test_merge_requires_identical_gene_universe(load_script_module):
    merge = load_script_module("merge_samples.py")

    merge.validate_gene_universe(["A", "B", "C"], ["C", "A", "B"], "sample_b")
    with pytest.raises(ValueError, match="1 missing and 1 extra"):
        merge.validate_gene_universe(
            ["A", "B", "C"],
            ["A", "B", "D"],
            "sample_b",
        )


def test_marker_filter_uses_direct_log2fc_threshold(monkeypatch, load_script_module):
    markers = load_script_module("07_markers_de.py")
    adata = clustered_adata()
    captured = {}

    monkeypatch.setattr(markers.sc.tl, "rank_genes_groups", lambda *args, **kwargs: None)

    def capture_filter(*args, **kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(markers.sc.tl, "filter_rank_genes_groups", capture_filter)
    markers.find_markers(
        adata,
        cluster_key="leiden",
        method="wilcoxon",
        n_genes=20,
        logfc_min=0.25,
        pct_min=0.1,
    )

    assert captured["min_fold_change"] == pytest.approx(0.25)
    assert adata.uns["marker_filter"]["min_log2_fold_change"] == pytest.approx(0.25)


def test_marker_scoring_preserves_unknown_for_small_margin(monkeypatch, load_script_module):
    annotation = load_script_module("08_annotate.py")
    adata = clustered_adata()

    def fake_score_genes(target, gene_list, score_name, **kwargs):
        target.obs[score_name] = (
            [1.0, 1.0, 2.0, 2.0]
            if gene_list[0].startswith("A")
            else [0.95, 0.95, 0.5, 0.5]
        )

    monkeypatch.setattr(annotation.sc.tl, "score_genes", fake_score_genes)
    _, summary, coverage = annotation.score_markers(
        adata,
        {"Type A": ["A1", "A2"], "Type B": ["B1", "B2"]},
        cluster_key="leiden",
        min_coverage=1.0,
        min_margin=0.1,
    )

    labels = summary.set_index("cluster")["marker_label"].to_dict()
    assert labels == {"0": "Unknown", "1": "Type A"}
    assert coverage["eligible"].all()
    assert set(adata.obs.loc[adata.obs["leiden"] == "0", "marker_cell_type"]) == {
        "Unknown"
    }


def test_manual_annotation_requires_every_cluster(tmp_path, load_script_module):
    annotation = load_script_module("08_annotate.py")
    mapping = tmp_path / "manual.tsv"
    pd.DataFrame({"cluster": ["0"], "cell_type": ["Type A"]}).to_csv(
        mapping, sep="\t", index=False
    )

    with pytest.raises(ValueError, match="missing clusters: 1"):
        annotation.load_manual(clustered_adata(), mapping, "leiden")


def test_manual_annotation_rejects_blank_labels(tmp_path, load_script_module):
    annotation = load_script_module("08_annotate.py")
    mapping = tmp_path / "manual.tsv"
    pd.DataFrame({"cluster": ["0", "1"], "cell_type": ["Type A", " "]}).to_csv(
        mapping, sep="\t", index=False
    )

    with pytest.raises(ValueError, match="cell_type cannot be blank"):
        annotation.load_manual(clustered_adata(), mapping, "leiden")


def test_custom_marker_dictionary_requires_two_unique_genes(tmp_path, load_script_module):
    annotation = load_script_module("08_annotate.py")
    markers = tmp_path / "markers.json"
    markers.write_text(json.dumps({"Type A": ["GENE1", "GENE1"]}), encoding="utf-8")

    with pytest.raises(ValueError, match="at least two unique genes"):
        annotation.load_marker_dictionary(markers)


def test_celltypist_input_is_rebuilt_from_raw_counts(load_script_module):
    annotation = load_script_module("08_annotate.py")
    adata = clustered_adata()
    adata.X = sparse.csr_matrix(np.full((4, 4), 99.0, dtype=np.float32))
    adata.layers["counts"] = sparse.csr_matrix(
        np.array(
            [
                [1, 2, 0, 1],
                [4, 0, 2, 2],
                [0, 3, 1, 0],
                [5, 1, 2, 2],
            ],
            dtype=np.float32,
        )
    )

    prepared = annotation.prepare_celltypist_input(adata)
    recovered_totals = np.asarray(np.expm1(prepared.X).sum(axis=1)).ravel()

    assert prepared.obs_names.equals(adata.obs_names)
    assert prepared.var_names.equals(adata.var_names)
    assert recovered_totals == pytest.approx(np.full(4, 10_000.0), rel=1e-5)
    assert not np.allclose(prepared.X.toarray(), adata.X.toarray())


def test_trajectory_barcode_map_is_one_to_one(tmp_path, load_script_module):
    trajectory = load_script_module("09_trajectory.py")
    loom = ad.AnnData(
        X=sparse.csr_matrix(np.ones((2, 3))),
        obs=pd.DataFrame(index=["loom_a", "loom_b"]),
        var=pd.DataFrame(index=["gene_a", "gene_b", "gene_c"]),
    )
    mapping = tmp_path / "barcode_map.tsv"
    pd.DataFrame(
        {
            "analysis_barcode": ["sample_a_cell", "sample_a_cell"],
            "loom_barcode": ["loom_a", "loom_b"],
        }
    ).to_csv(mapping, sep="\t", index=False)

    with pytest.raises(ValueError, match="one-to-one"):
        trajectory._apply_barcode_map(loom, mapping)


def test_trajectory_barcode_map_renames_matching_cells(tmp_path, load_script_module):
    trajectory = load_script_module("09_trajectory.py")
    loom = ad.AnnData(
        X=sparse.csr_matrix(np.ones((2, 3))),
        obs=pd.DataFrame(index=["loom_a", "loom_b"]),
        var=pd.DataFrame(index=["gene_a", "gene_b", "gene_c"]),
    )
    mapping = tmp_path / "barcode_map.tsv"
    pd.DataFrame(
        {
            "analysis_barcode": ["sample_a_cell", "sample_b_cell"],
            "loom_barcode": ["loom_a", "loom_b"],
        }
    ).to_csv(mapping, sep="\t", index=False)

    mapped = trajectory._apply_barcode_map(loom, mapping)
    assert mapped.obs_names.tolist() == ["sample_a_cell", "sample_b_cell"]


def test_trajectory_barcode_map_rejects_whitespace_only_ids(tmp_path, load_script_module):
    trajectory = load_script_module("09_trajectory.py")
    loom = ad.AnnData(
        X=sparse.csr_matrix(np.ones((2, 3))),
        obs=pd.DataFrame(index=["loom_a", "loom_b"]),
        var=pd.DataFrame(index=["gene_a", "gene_b", "gene_c"]),
    )
    mapping = tmp_path / "barcode_map.tsv"
    pd.DataFrame(
        {
            "analysis_barcode": ["sample_a_cell", "   "],
            "loom_barcode": ["loom_a", "loom_b"],
        }
    ).to_csv(mapping, sep="\t", index=False)

    with pytest.raises(ValueError, match="complete and one-to-one"):
        trajectory._apply_barcode_map(loom, mapping)


def test_trajectory_rejects_duplicate_loom_barcodes(load_script_module):
    trajectory = load_script_module("09_trajectory.py")
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        loom = ad.AnnData(
            X=sparse.csr_matrix(np.ones((2, 3))),
            obs=pd.DataFrame(index=["loom_a", "loom_a"]),
            var=pd.DataFrame(index=["gene_a", "gene_b", "gene_c"]),
        )

    with pytest.raises(ValueError, match="cell barcodes are duplicated"):
        trajectory._validate_loom_identifiers(loom)
