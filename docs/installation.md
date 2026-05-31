# Installation Guide

## Prerequisites

- Linux (Ubuntu 20.04+ recommended) or macOS
- 16 GB RAM minimum (32 GB+ for large datasets)
- ~10 GB disk space for environments

## 1. Install Miniforge (conda + mamba)

```bash
wget -O Miniforge3.sh \
  "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "${HOME}/miniforge3"
source "${HOME}/miniforge3/etc/profile.d/conda.sh"
conda init bash
source ~/.bashrc
```

## 2. Clone the repository

```bash
git clone https://github.com/your-org/scrna-pipeline.git
cd scrna-pipeline
```

## 3. Create conda environments

```bash
# Main Python environment (~5–8 min)
mamba env create -f envs/scrna.yaml

# Optional: R environment for scDblFinder / scran (~10–15 min)
mamba env create -f envs/r_env.yaml
```

## 4. Activate and verify

```bash
conda activate scrna

python -c "
import scanpy, anndata, scrublet, harmonypy
print('scanpy:', scanpy.__version__)
print('anndata:', anndata.__version__)
print('All imports OK')
"

snakemake --version
```

## 5. Run tests

```bash
pytest tests/ -v
```

## Troubleshooting

### "Package not found" during conda create
Update conda/mamba:
```bash
mamba update -n base mamba
```

### Slow environment creation
Use libmamba solver:
```bash
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

### rpy2 import error
Make sure you activated the correct environment:
```bash
conda activate scrna_r
python -c "import rpy2"
```

### Out of memory during pipeline
Reduce parallelism:
```bash
snakemake --use-conda --cores 2 --resources mem_mb=8000
```
