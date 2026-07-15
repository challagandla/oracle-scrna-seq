#!/usr/bin/env python3
"""
Merge per-sample AnnData objects into a single AnnData.
Adds sample / batch / condition metadata from samples.tsv.
"""

import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import pandas as pd
import scanpy as sc

from utils.io_utils import save_adata

sc.settings.verbosity = 1


def validate_gene_universe(reference_genes, observed_genes, sample):
    """Require the same named features in every sample; ordering may differ."""
    reference = set(map(str, reference_genes))
    observed = set(map(str, observed_genes))
    missing = sorted(reference - observed)
    extra = sorted(observed - reference)
    if missing or extra:
        missing_examples = ", ".join(missing[:5]) or "none"
        extra_examples = ", ".join(extra[:5]) or "none"
        raise ValueError(
            f"Sample {sample} has an incompatible gene universe relative to the "
            f"first sample: {len(missing)} missing and {len(extra)} extra genes. "
            f"Missing examples: {missing_examples}. Extra examples: {extra_examples}. "
            "All samples must use identical feature sets from the same reference "
            "annotation; gene order may differ."
        )


def main(args):
    os.makedirs(args.out, exist_ok=True)
    samples = pd.read_csv(args.samples_tsv, sep="\t", dtype=str)
    adatas = []
    reference_genes = None

    for _, row in samples.iterrows():
        sample = row["sample"]
        path   = os.path.join(args.input_dir, sample,
                              "01_doublets", "no_doublets.h5ad")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Required per-sample result is missing: {path}")
        adata = ad.read_h5ad(path)
        genes = adata.var_names.astype(str)
        if reference_genes is None:
            reference_genes = genes.copy()
        else:
            validate_gene_universe(reference_genes, genes, sample)
        # Propagate every sample-level design covariate. Input locations are not
        # biological metadata and are deliberately excluded from AnnData.
        metadata_columns = [
            column
            for column in samples.columns
            if column not in {"path", "raw_fastq_r1", "raw_fastq_r2"}
        ]
        for column in metadata_columns:
            value = row[column]
            if pd.notna(value):
                adata.obs[column] = str(value)
        # Prefix barcodes with sample name to avoid collisions
        adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
        adatas.append(adata)
        print(f"  Loaded {sample}: {adata.n_obs:,} cells")

    if not adatas:
        sys.exit("No samples loaded. Check paths and samples.tsv.")

    # Exact feature-set validation above makes this an order-safe concatenation
    # without inventing structural zeros for genes absent from a sample.
    merged = ad.concat(adatas, join="inner")
    merged.obs_names_make_unique()
    merged.var_names_make_unique()
    n_genes_before = merged.n_vars
    if args.min_cells > 0:
        sc.pp.filter_genes(merged, min_cells=args.min_cells)
    if merged.n_vars == 0:
        raise ValueError("Global gene filtering removed every gene")
    print(f"\nMerged: {merged.n_obs:,} cells × {merged.n_vars:,} genes")
    print(f"Global gene filter: {n_genes_before:,} → {merged.n_vars:,} genes")

    save_adata(merged, os.path.join(args.out, "merged.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--samples-tsv",  required=True)
    p.add_argument("--input-dir",    default="results")
    p.add_argument("--out",          default="results/merged/00_merged")
    p.add_argument("--min-cells",    type=int, default=20)
    main(p.parse_args())
