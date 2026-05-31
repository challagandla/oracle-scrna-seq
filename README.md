# 🧬 scRNA-seq Analysis Pipeline

[![CI](https://github.com/your-org/scrna-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/your-org/scrna-pipeline/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-green.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org)

End-to-end single-cell RNA-seq analysis pipeline following **[scverse best practices](https://www.sc-best-practices.org)**. Reproducible, containerised, and fully automated with Snakemake.

---

## 📋 Table of Contents

- [Pipeline overview](#pipeline-overview)
- [Quick start](#quick-start)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running the pipeline](#running-the-pipeline)
- [Step-by-step guide](#step-by-step-guide)
- [Output structure](#output-structure)
- [Customisation](#customisation)
- [Contributing](#contributing)
- [Citation](#citation)

---

## Pipeline Overview

```
Raw counts (.h5ad / 10X .h5 / MEX dir)
            │
  ┌─────────▼────────┐  per sample
  │  00_qc.py        │  MAD-based filtering, knee plot, scatter QC
  │  01_doublets.py  │  scrublet / scDblFinder
  └─────────┬────────┘
            │ merge_samples.py
            ▼
  ┌──────────────────┐  merged dataset
  │  02_normalize.py │  log-norm (or scran); raw counts preserved
  │  03_hvg_pca.py   │  HVG selection, PCA + elbow plot + loadings
  │  04_batch_*.py   │  Harmony / scVI / ComBat (skip if 1 batch)
  │  05_embedding.py │  kNN graph + UMAP
  │  06_cluster.py   │  Leiden @ multi-res + silhouette selection
  │  07_markers.py   │  Wilcoxon markers + pseudobulk DE (pydeseq2)
  │  08_annotate.py  │  Marker scoring + CellTypist + manual TSV
  │  09_trajectory.py│  PAGA + scVelo RNA velocity (optional)
  └──────────────────┘
            │
      report.html  (self-contained, all figures embedded)
```

---

## Quick Start

```bash
# 1. Clone repository
git clone https://github.com/your-org/scrna-pipeline.git
cd scrna-pipeline

# 2. Create conda environment
conda env create -f envs/scrna.yaml
conda activate scrna

# 3. Edit config and sample sheet
nano config/config.yaml
nano config/samples.tsv

# 4. Place raw data in data/
#    (or update paths in samples.tsv)

# 5. Run pipeline
snakemake --use-conda --cores 8

# 6. View report
open results/report.html   # or xdg-open on Linux
```

---

## Installation

### Requirements

| Tool | Version | Install |
|------|---------|---------|
| conda / mamba | any | [miniforge](https://github.com/conda-forge/miniforge) |
| Snakemake | ≥ 8.0 | included in `envs/scrna.yaml` |
| Python | 3.11 | included |
| R + scDblFinder | optional | `envs/r_env.yaml` |

### Step-by-step on Ubuntu

```bash
# Install Miniforge (skip if conda/mamba already installed)
wget -O Miniforge3.sh \
  "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "${HOME}/miniforge3"
source "${HOME}/miniforge3/etc/profile.d/conda.sh"
conda init bash && source ~/.bashrc

# Install mamba for faster solves
conda install -n base -c conda-forge mamba -y

# Create scRNA environment (takes ~5–10 min)
mamba env create -f envs/scrna.yaml

# (Optional) R environment for scDblFinder / scran
mamba env create -f envs/r_env.yaml

# Activate
conda activate scrna

# Verify
python -c "import scanpy; print(scanpy.__version__)"
snakemake --version
```

---

## Configuration

All analysis parameters live in **`config/config.yaml`**. The sample sheet lives in **`config/samples.tsv`**.

### Sample sheet (`config/samples.tsv`)

A tab-separated file with one row per sample:

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | ✅ | Unique sample identifier |
| `path` | ✅ | Path to `.h5ad`, `.h5`, or 10X MEX directory |
| `batch` | ✅ | Batch label (use same value if single batch) |
| `condition` | optional | Experimental condition for DE |
| `description` | optional | Free-text description |

Example:
```tsv
sample	path	batch	condition	description
pbmc_ctrl	data/pbmc_ctrl.h5ad	run1	control	PBMC donor 1
pbmc_treated	data/pbmc_treat.h5ad	run1	treated	PBMC donor 1, stimulated
```

### Key config parameters

```yaml
species: "human"          # or "mouse" — auto-sets MT/ribo/Hb prefixes

qc:
  mt_hard: 8              # raise to 10-15% for neurons/cardiomyocytes
  mad_counts: 5           # permissive; lower for high-quality datasets

batch:
  batch_key: "batch"      # set to "" to skip batch correction
  method: "harmony"       # harmony | scvi | combat

clustering:
  resolutions: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

annotation:
  celltypist: true        # automated cell type classification
```

Full parameter reference: [`docs/parameters.md`](docs/parameters.md)

---

## Running the Pipeline

### Full run (all steps)
```bash
snakemake --use-conda --cores 8
```

### Dry-run (show plan without executing)
```bash
snakemake -n --use-conda
```

### Run until a specific step
```bash
snakemake --use-conda --cores 8 --until cluster
```

### Run only specific samples
```bash
snakemake --use-conda --cores 4 results/pbmc_ctrl/00_qc/qc_filtered.h5ad
```

### Visualise DAG
```bash
snakemake --dag | dot -Tsvg > dag.svg
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

### HPC / cluster execution
```bash
# SLURM
snakemake --use-conda --cores 8 \
  --executor cluster-generic \
  --cluster-generic-submit-cmd "sbatch --cpus-per-task={threads} --mem={resources.mem_mb}M"

# SGE
snakemake --use-conda --cores 8 \
  --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem_mb}M"
```

### Override config on the command line
```bash
snakemake --use-conda --cores 8 \
  --config species=mouse \
           "batch={'batch_key':'batch','method':'scvi','n_latent':30}"
```

---

## Step-by-Step Guide

### Step 00 — Quality Control

**Purpose**: Remove low-quality cells (empty droplets, dying cells, debris) and lowly-expressed genes.

**Key outputs**:
- `figures/01_knee_plot.png` — Barcode-rank plot; the "knee" separates real cells from empty droplets
- `figures/02_qc_violin_before.png` — Distribution of all QC metrics pre-filtering
- `figures/03_qc_histograms.png` — Histograms with MAD-based threshold lines overlaid
- `figures/04_counts_vs_genes.png` — Scatter plot coloured by MT%; dying cells cluster bottom-left with high MT%
- `qc_stats.json` — Machine-readable before/after statistics

**Algorithm** (MAD = Median Absolute Deviation):
```
outlier_counts  = |log(counts) - median| > 5 × MAD  OR  counts < min_counts
outlier_genes   = |log(genes)  - median| > 5 × MAD  OR  genes  < min_genes
outlier_mt      = pct_mt > min(median + 3×MAD, mt_hard%)
fail_qc         = outlier_counts OR outlier_genes OR outlier_mt
```

> **Why MAD, not fixed thresholds?** Fixed cutoffs (e.g. "keep cells with >500 genes") fail because different tissues, protocols, and species produce different ranges. MAD adapts to your dataset. See [scverse best practices](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html).

**Adjust for your tissue**:
```yaml
qc:
  mt_hard: 15      # neurons, cardiomyocytes — naturally high MT%
  mad_mt: 5        # very permissive (studying rare stressed populations)
  min_counts: 200  # lower for SMART-seq / low-input data
```

---

### Step 01 — Doublet Detection

**Purpose**: Identify and remove doublets — droplets that captured 2+ cells. Doublets appear as fake "intermediate" cell types and inflate cluster counts.

**Methods**:

| Method | How | Accuracy | Requires |
|--------|-----|----------|----------|
| `scrublet` | Simulates doublets; computes similarity | Good | Python only |
| `scdblfinder` | Graph-based, uses known cell cycle | Best | R + Bioconductor |

**Recommended**: Use `scrublet` for speed; switch to `scdblfinder` for final analyses.

> **Rule of thumb**: 10X Genomics generates ~0.8% doublets per 1,000 cells captured. A 10,000-cell experiment expects ~8% doublets.

---

### Step 02 — Normalization

**Purpose**: Remove technical variation in library size (sequencing depth) so cells are comparable.

**Methods**:

| Method | Steps | Best for |
|--------|-------|----------|
| `log` (default) | Divide by library size → ×10,000 → log1p | Most datasets |
| `scran` | Pooling-based size factors | Low library size variation, many zeros |

Raw integer counts are always preserved in `adata.layers["counts"]` for:
- Pseudobulk differential expression (requires counts, not normalized values)
- scVI (deep generative model requires counts)
- Future re-analysis

`adata.raw` is frozen at the log-norm stage for marker gene testing.

---

### Step 03 — HVG Selection & PCA

**Purpose**: Select informative genes and compress expression space.

**HVG selection** (`seurat_v3` flavor):
- Uses raw counts to model mean-variance relationship
- Selects top 2,000 genes by normalized variance
- Batch-aware if `batch_key` is set (selects HVGs consistent across batches)
- Non-HVGs are retained in the object but excluded from PCA

**Scaling**: HVGs are zero-mean, unit-variance scaled (capped at 10) before PCA.
> Do NOT use scaled data for visualization or DE. Only for PCA.

**Choosing n_pcs** (from elbow plot):
- Look for where variance ratio stops decreasing sharply ("elbow")
- Rule of thumb: use PCs explaining 80% cumulative variance
- Typical range: 20–50 PCs

---

### Step 04 — Batch Correction

**Purpose**: Remove technical batch effects while preserving biological variation.

| Method | Pros | Cons |
|--------|------|------|
| Harmony | Fast (seconds), no retraining, works well | May over-integrate |
| scVI | Most powerful, probabilistic | Slow (minutes–hours), needs counts |
| ComBat | Simple, interpretable | Assumes linear effects only |

> **Skip if you have a single sample.** Over-integration is a real risk — if batches are already well-mixed, skip this step.

**Diagnostic**: Compare UMAP coloured by batch before vs after correction. After correction, cells should be interleaved without losing biological clusters.

---

### Step 05 — Neighbor Graph & UMAP

**Purpose**: Embed cells in 2D for visualization.

**kNN graph**: Finds `n_neighbors` nearest cells in PCA space. This graph is the foundation for both UMAP and Leiden clustering.

**UMAP parameters**:
- `n_neighbors` (15): higher → more global structure; lower → more local
- `min_dist` (0.3): lower → tighter clusters; higher → more spread

> ⚠️ **UMAP is for visualization only.** Distances in UMAP space are not biologically meaningful. Two clusters that look close may not be transcriptionally similar.

---

### Step 06 — Leiden Clustering

**Purpose**: Group cells into discrete clusters.

The pipeline runs Leiden at 6 resolutions and selects the best by **silhouette score** (measures how well-separated clusters are in PCA space).

> **Resolution guidance**:
> - Start with 0.4–0.8 for most datasets
> - More cells → can use higher resolution
> - Check: do clusters make biological sense? Is the silhouette curve clearly peaked?

Inspect `figures/01_umap_all_resolutions.png` to see how resolution affects cluster granularity.

---

### Step 07 — Marker Genes & DE

**A) Cluster markers** (Wilcoxon one-vs-rest):
- Non-parametric, robust to zero-inflation
- Recommended by scverse for scRNA-seq
- Results filtered by log2FC ≥ 0.25 and fraction expressing ≥ 10%

**B) Condition DE** (pseudobulk, requires biological replicates):
```yaml
markers:
  condition_key: "treatment"
  group1: "control"
  group2: "treated"
```
Pseudobulk aggregates counts per sample, then tests with DESeq2 — the only statistically valid approach when you have replicates.

> ❌ **Do not run per-cell DE between conditions** (inflated degrees of freedom). Always use pseudobulk.

---

### Step 08 — Cell Type Annotation

Three complementary approaches:

1. **Gene scoring** (`sc.tl.score_genes`): Scores each cell for each cell type using module genes. Fast, interpretable, works for any tissue with a known marker list.

2. **CellTypist**: ML classifier trained on 20+ human tissue atlases. Pass `celltypist: true` in config. Best for immune cells.

3. **Manual annotation**: After reviewing marker dotplots and UMAPs, create a TSV:
   ```tsv
   cluster	cell_type
   0	CD4 T cell
   1	B cell
   2	NK cell
   ```
   Then set `manual_tsv: "config/my_annotations.tsv"`.

**Composition plot**: `figures/06_composition_barplot.png` shows cell type proportions per sample — useful for identifying disease-associated compositional shifts.

---

### Step 09 — Trajectory (optional)

**PAGA** (any dataset):
```yaml
trajectory:
  run_paga: true
```
Computes cluster-level connectivity. A PAGA-initialized UMAP often preserves developmental trajectories better than random initialization.

**scVelo** (requires spliced/unspliced from velocyto or STARsolo):
```bash
# Generate .loom with velocyto (run after STAR/CellRanger alignment)
velocyto run-smartseq2 -o velocyto/ bam/*.bam genes.gtf
# or
velocyto run10x -m repeat_mask.gtf sample/ cellranger/refdata-gex/genes/genes.gtf
```
```yaml
trajectory:
  loom_file: "data/sample.loom"
```

---

## Output Structure

```
results/
├── report.html                          ← self-contained HTML report (open in browser)
│
├── <sample>/
│   ├── 00_qc/
│   │   ├── qc_filtered.h5ad
│   │   ├── qc_stats.json
│   │   └── figures/
│   │       ├── 01_knee_plot.png
│   │       ├── 02_qc_violin_before.png
│   │       ├── 03_qc_histograms.png      ← thresholds overlaid
│   │       ├── 04_counts_vs_genes.png
│   │       ├── 05_counts_vs_mt.png
│   │       └── 06_qc_violin_after.png
│   └── 01_doublets/
│       ├── no_doublets.h5ad
│       └── figures/
│           └── 01_doublet_score_hist.png
│
└── merged/
    ├── 00_merged/merged.h5ad
    ├── 02_normalize/
    │   ├── normalized.h5ad
    │   └── figures/{01,02}_count_distribution.png
    ├── 03_hvg_pca/
    │   ├── hvg_pca.h5ad
    │   └── figures/{01_hvg,02_elbow,03_pca_qc,04_loadings}.png
    ├── 04_batch/
    │   ├── batch_corrected.h5ad
    │   └── figures/{01_before,02_after}_correction.png
    ├── 05_embedding/
    │   ├── embedded.h5ad
    │   └── figures/{01_qc,02_batch,03_sample,04_condition}_umap.png
    ├── 06_cluster/
    │   ├── clustered.h5ad
    │   ├── clustering_summary.csv
    │   └── figures/{01_all_res,02_silhouette,03_selected,04_sizes}.png
    ├── 07_markers/
    │   ├── with_markers.h5ad
    │   ├── cluster_markers.csv
    │   ├── DE/pseudobulk_cluster_*.csv
    │   └── figures/{01_dotplot,02_heatmap,03_violin,04_volcano_*}.png
    ├── 08_annotate/
    │   ├── annotated.h5ad             ← final analysis-ready object
    │   └── figures/{01-06}_*.png
    └── 09_trajectory/
        ├── trajectory.h5ad
        ├── velocity.h5ad              (scVelo only)
        └── figures/{paga,velocity,latent_time}.png
```

---

## Customisation

### Custom marker genes (JSON)

```json
{
  "Microglia":      ["P2RY12", "TMEM119", "CX3CR1", "IBA1"],
  "Astrocyte":      ["GFAP", "AQP4", "ALDH1L1", "S100B"],
  "Oligodendrocyte":["MBP", "MOG", "PLP1", "MAG"],
  "Neuron":         ["RBFOX3", "MAP2", "TUBB3", "SYN1"]
}
```

Point to it in config:
```yaml
annotation:
  markers_json: "config/brain_markers.json"
```

### Mouse brain example config

```yaml
species: "mouse"
qc:
  mt_hard: 15       # neurons have higher MT%
  min_cells: 10     # for rare cell types
hvg_pca:
  n_top_genes: 3000
batch:
  batch_key: "sample"
  method: "harmony"
annotation:
  markers_json: "config/mouse_brain_markers.json"
trajectory:
  run_paga: true
```

---

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-improvement`
3. Commit your changes: `git commit -m "Add: my improvement"`
4. Push: `git push origin feature/my-improvement`
5. Open a Pull Request

Please run tests before submitting: `pytest tests/ -v`

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **scverse**: Virshup et al. (2023) *Nature Biotechnology* — anndata, scanpy
- **Leiden algorithm**: Traag et al. (2019) *Scientific Reports*
- **Harmony**: Korsunsky et al. (2019) *Nature Methods*
- **scrublet**: Wolock et al. (2019) *Cell Systems*
- **scDblFinder**: Germain et al. (2021) *F1000Research*
- **CellTypist**: Domínguez Conde et al. (2022) *Science*
- **scVelo**: Bergen et al. (2020) *Nature Biotechnology*
- **pydeseq2**: Muzellec et al. (2023) *Bioinformatics*
- **PAGA**: Wolf et al. (2019) *Genome Biology*
- **Best practices review**: Luecken & Theis (2019) *Molecular Systems Biology*

---

## License

MIT — see [LICENSE](LICENSE)
