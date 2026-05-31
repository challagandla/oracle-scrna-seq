#!/usr/bin/env python3
"""
Merge per-sample AnnData objects into a single AnnData.
Adds sample / batch / condition metadata from samples.tsv.
"""

import argparse, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import anndata as ad
import pandas as pd
import scanpy as sc

from utils.io_utils import save_adata

sc.settings.verbosity = 1


def main(args):
    os.makedirs(args.out, exist_ok=True)
    samples = pd.read_csv(args.samples_tsv, sep="\t", dtype=str)
    adatas = []

    for _, row in samples.iterrows():
        sample = row["sample"]
        path   = os.path.join(args.input_dir, sample,
                              "01_doublets", "no_doublets.h5ad")
        if not os.path.exists(path):
            print(f"WARNING: {path} not found — skipping {sample}")
            continue
        adata = ad.read_h5ad(path)
        # Tag each cell with sample metadata
        for col in ["sample", "batch", "condition"]:
            if col in row.index:
                adata.obs[col] = row[col]
        # Prefix barcodes with sample name to avoid collisions
        adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
        adatas.append(adata)
        print(f"  Loaded {sample}: {adata.n_obs:,} cells")

    if not adatas:
        sys.exit("No samples loaded. Check paths and samples.tsv.")

    merged = ad.concat(adatas, join="outer", fill_value=0)
    merged.obs_names_make_unique()
    merged.var_names_make_unique()
    print(f"\nMerged: {merged.n_obs:,} cells × {merged.n_vars:,} genes")

    save_adata(merged, os.path.join(args.out, "merged.h5ad"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--samples-tsv",  required=True)
    p.add_argument("--input-dir",    default="results")
    p.add_argument("--out",          default="results/merged/00_merged")
    main(p.parse_args())
