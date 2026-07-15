# ORACLE scRNA-seq matrix-to-report workflow

[![CI](https://github.com/challagandla/oracle-scrna-seq/actions/workflows/ci.yml/badge.svg)](https://github.com/challagandla/oracle-scrna-seq/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/Snakemake-8%2B-green.svg)](https://snakemake.readthedocs.io/)
[![Python](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)

ORACLE is a Snakemake workflow for reproducible downstream analysis of single-cell RNA-seq count matrices. It performs per-sample quality control and doublet detection, merges samples, normalizes counts, builds a low-dimensional representation, clusters cells at a configured biological resolution, ranks cluster markers, supports guarded pseudobulk differential expression, applies optional cell-type annotation, and produces a self-contained HTML report plus a run manifest.

This is deliberately a **count-matrix-to-report** workflow. Read alignment, UMI quantification, cell calling, demultiplexing, and ambient-RNA correction are upstream responsibilities. Optional FastQC, FastQ Screen, and MultiQC jobs form a parallel QC branch; they do not create or filter the expression matrix.

## Start here

The complete beginner-friendly guide is [TUTORIAL.md](TUTORIAL.md). It explains the biology and statistics as well as the commands, including experimental design, input validation, interpretation, troubleshooting, and limitations.

For a local run from the repository root:

```bash
conda env create -f envs/scrna.yaml
conda activate scrna

# Replace the example sample sheet and review every setting first.
$EDITOR config/samples.tsv
$EDITOR config/config.yaml

# Inspect the planned jobs before doing work.
snakemake --use-conda --cores 8 --dry-run

# Run and safely resume incomplete jobs if necessary.
snakemake --use-conda --cores 8 --rerun-incomplete --printshellcmds
```

When the workflow completes, open `results/report.html` and inspect `results/run_manifest.json`. A successful command is the start of review, not proof that the biological interpretation is correct.

## Input contract

Each row of `config/samples.tsv` represents one independent sample or capture. Every expression input must be one of:

- an `.h5ad` file whose `X` matrix contains nonnegative integer raw counts;
- a filtered 10x `.h5` count matrix; or
- a filtered 10x MEX directory.

All samples in one run must use the same species and the identical named gene
universe from one reference annotation; gene order may differ. The workflow
stops rather than turning absent features into false zero counts. Normalized,
log-transformed, scaled, or integrated values are not valid inputs. Unfiltered
droplet matrices must first undergo cell calling outside this workflow.

For condition-level inference, `sample` must identify the biological replicate. Cells are observations within a replicate; they are not independent replicates.

## Workflow

```text
validated raw count matrices
        |
        +-- optional raw-FASTQ QC branch (independent; no quantification)
        |
        v
per-sample cell QC -> per-sample doublet detection
        |
        v
merge -> normalize -> HVGs/PCA -> optional batch correction
        |
        v
neighbors/UMAP -> Leiden diagnostics + configured resolution
        |
        v
cluster markers -> guarded pseudobulk DE -> optional annotation
        |
        +-- optional PAGA / scVelo
        |
        v
HTML report + tables + H5AD objects + run manifest
```

## Scientific guardrails

- QC thresholds are starting points, not universal biological truths. Review every sample's distributions and retention.
- Correct only genuine technical batches. If batch and condition are confounded, integration cannot recover the missing design.
- UMAP is a visualization; its distances, islands, and apparent paths are not statistical evidence.
- The configured Leiden resolution is the reported clustering. Silhouette values across candidate resolutions are diagnostics, not an automatic definition of cell types.
- Cluster-marker tests rank features for annotation. Their cell-level p-values are not condition-level inference.
- Condition DE uses raw-count pseudobulks at the biological-sample level. Completed runs record per-cluster outcomes in `pseudobulk_status.csv`; if no requested comparison is eligible, the rule fails and prints the status table to its durable log.
- Annotation is opt-in. An incompatible marker dictionary fails fast at the coverage check; weak cluster scores, score margins, or classifier confidence remain `Unknown`. A manual map must cover every cluster and may explicitly assign `Unknown`.
- Cell-type composition tables are descriptive. Replicate-aware differential-abundance testing is outside the workflow.
- PAGA and RNA velocity are optional hypothesis-generating analyses with strict biological and identifier requirements.

## Main outputs at the default paths

| Output | Meaning |
|---|---|
| `results/report.html` | Self-contained visual summary for review |
| `results/run_config.yaml` | Effective configuration after Snakemake overrides |
| `results/run_manifest.json` | Configuration/input fingerprints, Conda-spec hashes, detected report-environment packages, Git SHA, seed, and warnings |
| `results/<sample>/00_qc/` | Per-sample QC object, statistics, and figures |
| `results/<sample>/01_doublets/` | Per-sample singlet object and doublet diagnostic |
| `results/merged/06_cluster/` | Clustered object and resolution diagnostics |
| `results/merged/07_markers/` | Cluster markers, pseudobulk results, and DE status table |
| `results/merged/08_annotate/` | Annotated object, annotation summary, and composition table |
| `results/merged/09_trajectory/` | Optional PAGA and velocity artifacts |
| `results/multiqc/multiqc_report.html` | Optional raw-FASTQ QC summary or explicit no-FASTQ notice |
| `logs/` | One log per rule; inspect these first when a job fails or is skipped |

## Reference documentation

- [Complete tutorial](TUTORIAL.md)
- [Single-cell best-practices reference](https://www.sc-best-practices.org/)

## License

MIT — see [LICENSE](LICENSE).
