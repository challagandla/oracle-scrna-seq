#!/usr/bin/env python3
"""Create deterministic raw-count inputs for the workflow smoke test."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def main(output_directory):
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260714)
    genes = ["MT-ND1", "MT-CO1", "RPL3", "RPS6"] + [f"GENE{index:04d}" for index in range(396)]
    rows = []
    design = [
        ("control_1", "control", "donor1", "batch1"),
        ("control_2", "control", "donor2", "batch2"),
        ("treated_1", "treated", "donor1", "batch1"),
        ("treated_2", "treated", "donor2", "batch2"),
    ]
    for sample, condition, donor, batch in design:
        counts = rng.negative_binomial(5, 0.5, size=(100, len(genes))).astype(np.int32)
        if condition == "treated":
            counts[:, 20:30] += rng.poisson(2, size=(100, 10)).astype(np.int32)
        adata = ad.AnnData(
            X=sparse.csr_matrix(counts),
            obs=pd.DataFrame(index=[f"{sample}_cell_{index:03d}" for index in range(100)]),
            var=pd.DataFrame(index=genes),
        )
        path = output / f"{sample}.h5ad"
        adata.write_h5ad(path)
        rows.append(
            {
                "sample": sample,
                "path": str(path.resolve()),
                "species": "human",
                "batch": batch,
                "condition": condition,
                "donor": donor,
                "description": "deterministic synthetic smoke-test input",
                "raw_fastq_r1": "",
                "raw_fastq_r2": "",
            }
        )
    pd.DataFrame(rows).to_csv(output / "samples.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=".test/data")
    main(parser.parse_args().out)
