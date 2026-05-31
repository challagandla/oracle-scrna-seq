# Full Parameter Reference

All parameters are set in `config/config.yaml`.

## Top-level

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samples_tsv` | `config/samples.tsv` | Path to sample sheet |
| `species` | `human` | `human` or `mouse`; auto-sets gene prefixes |
| `mt_prefix` | `""` | Override MT prefix (`""` = auto from species) |
| `ribo_prefix` | `""` | Override ribosomal prefix |
| `hb_pattern` | `""` | Override Hb gene regex |

## qc

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mad_counts` | `5` | MADs for total counts (log-scale). Lower = stricter. |
| `mad_genes` | `5` | MADs for gene counts (log-scale). |
| `mad_mt` | `3` | MADs for MT% (raw). More stringent than counts. |
| `mt_hard` | `8` | Hard MT% ceiling. Raise to 10–15 for neurons. |
| `min_cells` | `20` | Gene must appear in ≥ N cells. |
| `min_counts` | `500` | Hard floor for total counts per cell. |
| `min_genes` | `200` | Hard floor for detected genes per cell. |

## doublets

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `scrublet` | `scrublet` or `scdblfinder` |
| `expected_rate` | `0.06` | Expected doublet fraction per lane (10X: ~0.8% per 1k cells). |

## normalization

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `log` | `log` (library-size + log1p) or `scran` (pooling) |
| `target_sum` | `10000` | Target counts/cell for log-norm. |

## hvg_pca

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_top_genes` | `2000` | HVGs to select. 1500–3000 is typical. |
| `flavor` | `seurat_v3` | HVG method: `seurat_v3`, `seurat`, `cell_ranger`. |
| `n_pcs` | `50` | PCs to compute. |
| `n_pcs_use` | `30` | PCs to use downstream. See elbow plot. |

## batch

| Parameter | Default | Description |
|-----------|---------|-------------|
| `batch_key` | `""` | obs column with batch labels. `""` = skip. |
| `method` | `harmony` | `harmony`, `scvi`, or `combat` |
| `n_latent` | `30` | scVI latent dimensions. |

## embedding

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_neighbors` | `15` | kNN graph neighbors. 10–50. |
| `metric` | `euclidean` | Distance metric for kNN. |
| `min_dist` | `0.3` | UMAP: controls cluster tightness. 0.1–0.8. |
| `spread` | `1.0` | UMAP: scale of embedded space. |
| `use_rep` | `""` | Auto-set from batch method. Override if needed. |

## clustering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `resolutions` | `[0.2, 0.4, 0.6, 0.8, 1.0, 1.2]` | Leiden resolutions to test. |
| `default_resolution` | `0.6` | Fallback if silhouette fails. |

## markers

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `wilcoxon` | `wilcoxon`, `t-test`, `logreg` |
| `n_genes` | `50` | Top N marker genes to store per cluster. |
| `logfc_min` | `0.25` | Minimum log2 fold-change for marker filtering. |
| `pct_min` | `0.1` | Minimum fraction of cells expressing gene. |
| `condition_key` | `""` | obs column for condition DE. `""` = skip. |
| `group1` | `""` | Reference group label. |
| `group2` | `""` | Test group label. |

## annotation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `markers_json` | `""` | Path to `{cell_type: [gene1,...]}` JSON. `""` = built-in PBMC. |
| `celltypist` | `false` | Run CellTypist automated classification. |
| `celltypist_model` | `Immune_All_Low.pkl` | CellTypist model name. See [models](https://www.celltypist.org/models). |
| `manual_tsv` | `""` | Path to cluster→cell_type TSV. `""` = skip. |

## trajectory

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_paga` | `false` | Run PAGA cluster-level trajectory. |
| `loom_file` | `""` | Path to velocyto .loom for scVelo. `""` = skip. |
