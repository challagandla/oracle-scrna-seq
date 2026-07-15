"""Validation helpers for biological matrix inputs."""

from __future__ import annotations

import numpy as np
from scipy import sparse


def validate_raw_counts(adata, *, context: str = "input", atol: float = 1e-6) -> None:
    """Fail fast unless ``adata.X`` is a usable raw, non-negative count matrix."""
    if adata.n_obs < 3:
        raise ValueError(f"{context}: at least 3 cells are required; found {adata.n_obs}.")
    if adata.n_vars < 3:
        raise ValueError(f"{context}: at least 3 genes are required; found {adata.n_vars}.")
    if not adata.obs_names.is_unique:
        raise ValueError(f"{context}: cell barcodes are not unique.")
    if not adata.var_names.is_unique:
        raise ValueError(
            f"{context}: gene names are not unique. Make names unique before analysis."
        )

    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X).ravel()
    if values.size == 0 or not np.any(values != 0):
        raise ValueError(f"{context}: the count matrix contains no non-zero values.")
    if not np.isfinite(values).all():
        raise ValueError(f"{context}: the count matrix contains NaN or infinite values.")
    if np.min(values) < 0:
        raise ValueError(f"{context}: raw counts cannot be negative.")
    if not np.allclose(values, np.rint(values), atol=atol, rtol=0):
        raise ValueError(
            f"{context}: adata.X must contain raw integer counts, not normalized or "
            "log-transformed expression."
        )
