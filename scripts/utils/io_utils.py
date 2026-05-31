"""I/O helpers — loading .h5ad and 10X .h5 files."""
import os
import anndata as ad
import scanpy as sc


def load_adata(path: str) -> ad.AnnData:
    """Auto-detect format and load AnnData."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".h5ad":
        adata = ad.read_h5ad(path)
    elif ext == ".h5":
        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
    elif os.path.isdir(path):
        # 10X MEX directory (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
        adata.var_names_make_unique()
    else:
        raise ValueError(f"Unsupported input format: {ext}. Use .h5ad, .h5, or MEX dir.")
    print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes  [{path}]")
    return adata


def save_adata(adata: ad.AnnData, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    adata.write_h5ad(path)
    print(f"Saved:  {path}")
