"""Regression tests for project preflight and raw-count input contracts."""

from __future__ import annotations

from copy import deepcopy

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from utils.validation import validate_raw_counts
from validate_config import validate_project


def valid_config():
    return {
        "species": "human",
        "raw_fastq_qc": {
            "enabled": False,
            "fastq_screen_reads": ["R2"],
        },
        "qc": {"min_cells": 0, "min_counts": 0, "min_genes": 0},
        "hvg_pca": {"n_pcs": 20, "n_pcs_use": 10},
        "batch": {"batch_key": "", "method": "none"},
        "embedding": {"use_rep": ""},
        "clustering": {"resolutions": [0.5], "default_resolution": 0.5},
        "markers": {
            "method": "wilcoxon",
            "condition_key": "",
            "group1": "",
            "group2": "",
            "covariates": [],
            "min_replicates_per_group": 2,
        },
        "annotation": {"marker_set": "", "markers_json": "", "manual_tsv": ""},
        "trajectory": {"loom_file": "", "barcode_map_tsv": ""},
    }


def valid_samples(tmp_path):
    paths = []
    for sample in ("sample_a", "sample_b"):
        path = tmp_path / f"{sample}.h5ad"
        path.touch()
        paths.append(str(path))
    return pd.DataFrame(
        {
            "sample": ["sample_a", "sample_b"],
            "path": paths,
            "species": ["human", "human"],
            "batch": ["batch_1", "batch_2"],
            "condition": ["control", "treated"],
        }
    )


def test_project_preflight_accepts_valid_minimal_project(tmp_path):
    validate_project(valid_config(), valid_samples(tmp_path))


@pytest.mark.parametrize("identifier", ["../escape", "has whitespace", "/absolute"])
def test_project_preflight_rejects_unsafe_sample_ids(tmp_path, identifier):
    samples = valid_samples(tmp_path)
    samples.loc[0, "sample"] = identifier
    with pytest.raises(ValueError, match="sample IDs must match"):
        validate_project(valid_config(), samples)


def test_project_preflight_rejects_duplicate_sample_ids(tmp_path):
    samples = valid_samples(tmp_path)
    samples.loc[1, "sample"] = samples.loc[0, "sample"]
    with pytest.raises(ValueError, match="sample IDs must be unique"):
        validate_project(valid_config(), samples)


def test_project_preflight_rejects_empty_sample_sheet(tmp_path):
    samples = valid_samples(tmp_path).iloc[0:0]
    with pytest.raises(ValueError, match="at least one sample"):
        validate_project(valid_config(), samples)


def test_project_preflight_rejects_partial_fastq_pair(tmp_path):
    samples = valid_samples(tmp_path)
    read = tmp_path / "sample_a_R1.fastq.gz"
    read.touch()
    samples["raw_fastq_r1"] = [str(read), ""]
    samples["raw_fastq_r2"] = ["", ""]
    with pytest.raises(ValueError, match="only one raw FASTQ mate"):
        validate_project(valid_config(), samples)


def test_project_preflight_rejects_mixed_species(tmp_path):
    samples = valid_samples(tmp_path)
    samples.loc[1, "species"] = "mouse"
    with pytest.raises(ValueError, match="all samples must match global species"):
        validate_project(valid_config(), samples)


def test_project_preflight_enforces_pseudobulk_replicate_floor(tmp_path):
    config = deepcopy(valid_config())
    config["markers"].update(
        {
            "condition_key": "condition",
            "group1": "control",
            "group2": "treated",
            "min_replicates_per_group": 2,
        }
    )
    with pytest.raises(ValueError, match="fewer than 2 samples"):
        validate_project(config, valid_samples(tmp_path))


def test_project_preflight_rejects_batch_condition_confounding(tmp_path):
    config = deepcopy(valid_config())
    config["batch"].update({"batch_key": "batch", "method": "harmony"})

    with pytest.raises(ValueError, match="perfectly confounded with 'condition'"):
        validate_project(config, valid_samples(tmp_path))


def test_project_preflight_rejects_incomplete_barcode_map(tmp_path):
    config = deepcopy(valid_config())
    loom = tmp_path / "velocity.loom"
    loom.touch()
    barcode_map = tmp_path / "barcode_map.tsv"
    pd.DataFrame(
        {
            "analysis_barcode": ["sample_a_cell_1", ""],
            "loom_barcode": ["cell_1", "cell_2"],
        }
    ).to_csv(barcode_map, sep="\t", index=False)
    config["trajectory"].update(
        {"loom_file": str(loom), "barcode_map_tsv": str(barcode_map)}
    )
    with pytest.raises(ValueError, match="complete and one-to-one"):
        validate_project(config, valid_samples(tmp_path))


def make_adata(matrix, obs_names=None):
    matrix = np.asarray(matrix)
    return ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            index=obs_names or [f"cell_{index}" for index in range(matrix.shape[0])]
        ),
        var=pd.DataFrame(index=[f"gene_{index}" for index in range(matrix.shape[1])]),
    )


def test_raw_count_validation_accepts_integer_sparse_matrix():
    adata = make_adata(np.arange(12).reshape(4, 3))
    adata.X = sparse.csr_matrix(adata.X)
    validate_raw_counts(adata)


@pytest.mark.parametrize(
    ("matrix", "message"),
    [
        (np.full((4, 3), 0.5), "raw integer counts"),
        (np.full((4, 3), -1), "cannot be negative"),
        (np.zeros((4, 3)), "contains no non-zero values"),
    ],
)
def test_raw_count_validation_rejects_invalid_dense_matrices(matrix, message):
    with pytest.raises(ValueError, match=message):
        validate_raw_counts(make_adata(matrix))


def test_raw_count_validation_rejects_duplicate_barcodes():
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        adata = make_adata(
            np.ones((4, 3)), obs_names=["cell", "cell", "cell_2", "cell_3"]
        )
    with pytest.raises(ValueError, match="barcodes are not unique"):
        validate_raw_counts(adata)
