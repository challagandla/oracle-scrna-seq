# ORACLE scRNA-seq: a beginner-to-advanced matrix-to-report tutorial

This tutorial takes you from the biological question and count-matrix inputs to a reviewed HTML report. It is written for a first-time user, but it does not hide the statistical decisions that determine whether a result is trustworthy.

The central rule is simple:

> A workflow can complete without the experiment supporting the conclusion. Review the study design, diagnostics, warnings, and effect sizes before interpreting any result.

## Contents

1. [What single-cell RNA-seq measures](#1-what-single-cell-rna-seq-measures)
2. [What this workflow does](#2-what-this-workflow-does)
3. [Design the experiment before the analysis](#3-design-the-experiment-before-the-analysis)
4. [Install the software](#4-install-the-software)
5. [Prepare valid count-matrix inputs](#5-prepare-valid-count-matrix-inputs)
6. [Build the sample sheet](#6-build-the-sample-sheet)
7. [Configure the analysis](#7-configure-the-analysis)
8. [Plan and run the workflow](#8-plan-and-run-the-workflow)
9. [Interpret every analysis stage](#9-interpret-every-analysis-stage)
10. [Understand the output files](#10-understand-the-output-files)
11. [Troubleshoot safely](#11-troubleshoot-safely)
12. [Know the workflow's limits](#12-know-the-workflows-limits)
13. [Final review checklist](#13-final-review-checklist)
14. [Glossary](#14-glossary)

## 1. What single-cell RNA-seq measures

A single-cell RNA-seq expression matrix has cells as rows and genes as columns. Each value is usually a nonnegative integer count of unique molecules assigned to one cell and one gene.

The matrix contains several kinds of variation at once:

- biological cell identity, state, activation, and differentiation;
- differences among donors or other biological replicates;
- library size and transcript-capture efficiency;
- damaged cells, empty droplets retained by upstream cell calling, and doublets;
- sequencing run, chemistry, handling, and other technical batches;
- ambient RNA, alignment, reference-annotation, and quantification effects inherited from upstream processing.

Downstream analysis tries to preserve the biological signal while identifying or controlling technical variation. It cannot repair a confounded experiment or recover information that was never measured.

### The questions this workflow can help answer

The workflow can help you explore questions such as:

- Which transcriptional populations are present?
- Which genes characterize a cluster?
- Are similar cell states distributed across samples and batches?
- Within an annotated population, are there replicate-supported expression differences between two conditions?
- What descriptive differences in cell-type composition are visible among samples?
- Is a continuous process plausible enough to justify an optional trajectory analysis?

The wording matters. A UMAP island is not automatically a cell type. A marker p-value is not automatically condition-level evidence. A difference in cell percentages is not automatically differential abundance.

## 2. What this workflow does

ORACLE is a **count-matrix-to-report** workflow. Its main path is:

```text
raw integer count matrices
        |
        v
per-sample QC -> per-sample doublet detection
        |
        v
merge -> normalize -> select HVGs -> PCA
        |
        v
optional batch correction -> neighbors -> UMAP
        |
        v
Leiden resolution diagnostics -> configured clustering
        |
        v
cluster markers -> guarded pseudobulk DE
        |
        v
optional annotation -> optional trajectory
        |
        v
HTML report + H5AD objects + tables + run manifest
```

### Optional raw-FASTQ QC is a parallel branch

When paired FASTQ paths are supplied, Snakemake can also run FastQC, FastQ Screen, and MultiQC. This branch is scheduled alongside the count-matrix analysis. It does not:

- align reads;
- quantify genes or UMIs;
- call cell-containing droplets;
- create the matrix used by the main workflow; or
- automatically accept or reject cells based on a FastQC result.

For common 10x gene-expression libraries, the cDNA read is normally the informative read for contamination screening. The default `raw_fastq_qc.fastq_screen_reads: [R2]` reflects that convention, but you must confirm the read structure for your chemistry.

### Upstream responsibilities

Complete these before using the workflow:

- demultiplexing and sample assignment;
- chemistry-aware alignment and UMI quantification;
- cell calling from an unfiltered droplet matrix;
- ambient-RNA assessment or correction, when appropriate;
- reference-genome and gene-annotation selection;
- donor demultiplexing or genotype-based identity checks, when needed.

Keep the upstream reference, software versions, commands, and cell-calling decisions with the study provenance. The downstream run manifest cannot reconstruct upstream choices that it never sees.

## 3. Design the experiment before the analysis

No downstream method can separate two effects that were never separated experimentally.

### Cell, capture, sample, and biological replicate

These terms are not interchangeable:

| Term | Meaning |
|---|---|
| Cell | One measured cell or nucleus; cells from the same sample are not independent biological replicates |
| Capture/library | One library-preparation unit, such as one 10x channel |
| Sample | The unit named in `config/samples.tsv`; for pseudobulk inference it must represent an independent biological replicate |
| Biological replicate | An independently sampled organism, donor, culture, tissue unit, or other experimental unit appropriate to the question |
| Technical replicate | Repeated measurement of the same biological material; it does not increase biological replication |
| Condition | The biological group to compare, such as control versus treated |
| Batch | A technical grouping, such as processing day or sequencing run |

Thousands of cells from one donor still provide one donor for donor-level inference.

### Replication

The pseudobulk guardrail uses `markers.min_replicates_per_group` to refuse unsupported tests. Passing that minimum is not a power calculation and does not guarantee a reliable study. More independent replicates usually improve variance estimation, robustness to outliers, and the ability to model covariates.

Before collecting data, decide:

- what the experimental unit is;
- how many independent units are in each condition;
- whether measurements are paired or repeated;
- which covariates are scientifically required;
- how samples will be balanced across technical batches.

### Confounding

A design is confounded when two effects always occur together. For example:

```text
bad design:  every control is processed in batch 1
             every treated sample is processed in batch 2
```

In that design, a condition difference cannot be distinguished from a batch difference. Harmony, scVI, ComBat, or any other correction cannot recreate the missing comparison.

A better design distributes each condition across batches:

```text
batch 1: controls and treated samples
batch 2: controls and treated samples
```

Use a batch key only for a genuine technical effect. Do not label donor, condition, tissue, disease, or developmental stage as “batch” merely because those groups separate in PCA or UMAP.

### Paired designs and covariates

Optional `markers.covariates` names sample-sheet columns to include in the pseudobulk design. A paired study might use a donor column as a covariate. Every covariate must vary in a way that leaves the design estimable; a covariate that perfectly predicts condition is confounded.

The workflow supports guarded two-group pseudobulk analysis. Complex longitudinal, nested, random-effect, or multi-factor designs require a statistician and a method designed for that structure.

## 4. Install the software

### Practical requirements

- Linux is the primary execution platform.
- Git is required to clone the repository.
- Conda or Mamba is required for isolated rule environments.
- Python 3.11 and Snakemake 8 or later are defined by the environment files.
- 16 GB RAM is a practical minimum for small datasets; 32 GB or more is safer for multi-sample studies.
- Allow substantial disk space. H5AD intermediates preserve counts and analysis state and can be much larger than the original matrix.

A GPU is optional. It can accelerate selected scVI workloads but is not required for the standard PCA/Harmony path.

### Create the command environment

If you do not already have the repository, obtain it and enter its root directory:

```bash
git clone https://github.com/challagandla/oracle-scrna-seq.git
cd oracle-scrna-seq
```

Then create and activate the command environment:

```bash
conda env create -f envs/scrna.yaml
conda activate scrna

python -c "import scanpy; print('scanpy', scanpy.__version__)"
snakemake --version
```

Mamba may replace `conda` in the creation command if it is available.

Snakemake's `--use-conda` option creates or reuses the environment attached to each rule. The optional scDblFinder and scran branches use `envs/r_env.yaml`; raw-FASTQ QC uses `envs/raw_fastq_qc.yaml`.

### What the automated smoke test covers

The repository's deterministic tiny workflow test exercises the default Scrublet, log-normalization, no-batch-correction, no-annotation, and no-trajectory path. It checks workflow wiring and generic output creation; synthetic data cannot validate a biological conclusion.

The declared scDblFinder, scran, Harmony, scVI, ComBat, CellTypist, raw-FASTQ QC, PAGA, and scVelo branches require branch-specific validation with compatible inputs before study use. Availability in the configuration is not evidence that a particular reference, chemistry, dataset size, or experimental design is appropriate.

### Reproducibility boundary

The YAML files make environments isolated and recreatable from declared constraints, but broad version constraints are not a bit-for-bit lockfile. The manifest hashes every declared Conda YAML and records package versions detectable in the report rule's environment; it does not enumerate the resolved R or raw-QC environments. For a study that needs exact reconstruction, also export the fully solved environments or preserve explicit lockfiles.

## 5. Prepare valid count-matrix inputs

### Supported formats

Each sample path must point to one of:

- `.h5ad`: an AnnData object with cells in rows, genes in columns, and raw counts in `X`;
- `.h5`: a filtered 10x gene-expression HDF5 matrix; or
- a filtered 10x MEX directory containing its matrix, barcodes, and features files.

Use filtered, cell-called matrices. An unfiltered droplet matrix contains many background barcodes and requires an upstream cell-calling method.

### Raw integer counts are mandatory

The input `X` matrix must contain finite, nonnegative integer molecule counts. Do not supply:

- log-normalized values;
- counts per million or transcripts per million;
- scaled or z-scored expression;
- batch-corrected expression;
- an integrated embedding; or
- a matrix reconstructed by rounding normalized values.

The workflow needs real counts for QC, doublet detection, Seurat-v3 HVG selection, scVI, and pseudobulk DE. It preserves those counts in `layers["counts"]` before changing `X` during normalization.

If a downloaded `.h5ad` stores counts in a layer while `X` is normalized, create a new input whose `X` contains the original counts. Do not silently substitute a layer without recording what you changed.

### Gene identifiers

Gene names must be unique. Every sample in one run must contain the identical
set of named genes from the same reference annotation, although feature order may
differ. The workflow checks this at merge time and stops with examples of
missing and extra genes. It never fills genuinely absent features with zero,
because those structural zeros can become false biological signal. Gene symbols
are needed for the built-in marker sets and species-specific mitochondrial,
ribosomal, and hemoglobin patterns.

Do not mix, for example:

- human symbols in one sample and Ensembl IDs in another;
- gene symbols from different annotation releases without reconciliation;
- human and mouse matrices in one run; or
- gene IDs and protein names in a marker file.

If you use Ensembl IDs, explicitly map them to stable, unique symbols or provide compatible marker and QC definitions. Record the mapping version.

### Species consistency

One run must contain one species. The `species` setting controls gene-prefix defaults; it does not translate orthologs or make mixed-species integration valid.

### Cell barcodes

Within each sample, cell IDs must be unique. The merge stage prefixes barcodes with the sample ID, yielding analysis IDs such as:

```text
sample01_AAACCCAAGAAACACT-1
```

This convention matters later if RNA velocity is enabled.

### Preflight behavior

The workflow validates the sample sheet and inputs before expensive analysis. Treat a preflight failure as a data-contract problem. Do not bypass it by rounding counts, renaming species, or weakening checks without understanding the cause.

## 6. Build the sample sheet

Edit `config/samples.tsv` as a true tab-separated file. Do not align columns with spaces.

### Core columns

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Unique, filesystem-safe sample identifier and pseudobulk replicate ID |
| `path` | yes | Count-matrix path relative to the repository root or an absolute path |
| `species` | yes | `human` or `mouse`; every row must match the global `species` setting |
| `batch` | yes | Technical batch label; use one shared value if there is no batch effect to model |
| `condition` | yes | Biological comparison group; use one shared label when condition DE is disabled |
| `description` | no | Human-readable note |
| `raw_fastq_r1` | no | Optional read-1 FASTQ path |
| `raw_fastq_r2` | no | Optional read-2 FASTQ path |

Additional columns may define sample-level covariates used by `markers.covariates`.

Sample identifiers must match `[A-Za-z0-9][A-Za-z0-9_.-]*`: begin with a letter or number, then use only letters, numbers, dots, underscores, or hyphens. Whitespace, slashes, shell characters, and duplicate names are rejected.

### A coherent example

This fictional example contains one species, independent donors, and both conditions in both technical batches:

```tsv
sample	path	species	batch	condition	description
donor01	data/donor01_filtered.h5	human	day1	control	independent donor 01
donor02	data/donor02_filtered.h5	human	day2	control	independent donor 02
donor03	data/donor03_filtered.h5	human	day1	treated	independent donor 03
donor04	data/donor04_filtered.h5	human	day2	treated	independent donor 04
```

The minimum replication guard is not a claim that two samples per group are sufficient for every biological question. Plan replication using expected variability, effect size, tissue, cell abundance, and the intended statistical model.

### Technical lanes versus biological samples

If several sequencing lanes measure the same library, do not present those lanes as independent `sample` values for DE. Combine technical data upstream or use a documented strategy that preserves the biological replicate as the inferential unit.

### Raw FASTQ pairs

Supply both FASTQ columns for a sample or leave both blank. Confirm whether R1 or R2 contains cDNA for the assay. `raw_fastq_qc.fastq_screen_reads` controls which reads are screened; the default `[R2]` is appropriate for many, but not all, 10x gene-expression chemistries.

## 7. Configure the analysis

Edit `config/config.yaml` deliberately. The shipped values are starting points, not recommendations for every tissue or protocol.

### Paths, seed, and provenance

```yaml
samples_tsv: "config/samples.tsv"
results_dir: "results"
logs_dir: "logs"
random_seed: 42

provenance:
  hash_inputs: false
```

This tutorial shows paths under the default `results` and `logs` directories. If you change either root, substitute the configured location everywhere below.

The seed controls supported stochastic workflow steps, but a seed does not guarantee bit-for-bit agreement across different operating systems, hardware, thread counts, or library versions.

The run manifest always records absolute input paths, size, and modification time, plus SHA-256 hashes of the effective configuration, sample sheet, and every declared Conda YAML. Set `provenance.hash_inputs: true` when the manifest must also hash the full contents of every input matrix and FASTQ. Content hashing is stronger provenance but can add substantial I/O time for large studies.

### Optional raw-FASTQ QC

```yaml
raw_fastq_qc:
  enabled: true
  fastq_screen_conf: "config/fastq_screen.conf"
  fastq_screen_subset: 100000
  fastq_screen_reads: [R2]
```

If FASTQs are present, copy `config/fastq_screen.conf.example` to the configured local file and replace every database placeholder with a valid Bowtie2 index prefix. A FastQ Screen database path is site-specific and should not be treated as a portable reference definition.

When no FASTQs are supplied, the branch is skipped and the MultiQC page states that no raw reads were configured. This is not evidence that raw-read QC passed.

### Species and QC gene patterns

```yaml
species: "human"  # or "mouse"
mt_prefix: ""
ribo_prefix: ""
hb_pattern: ""
```

Blank patterns use species defaults. After the first QC run, verify that mitochondrial, ribosomal, and hemoglobin genes were actually detected. Zero detected mitochondrial genes often indicates incompatible gene identifiers, not exceptionally healthy cells.

### Cell QC

```yaml
qc:
  mad_counts: 5
  mad_genes: 5
  mad_mt: 3
  mt_hard: 8
  min_cells: 20
  min_counts: 500
  min_genes: 200
```

The workflow combines data-adaptive median-absolute-deviation bounds with hard floors or ceilings. Review the distributions before changing values.

Do not choose a threshold only because it retains a preferred percentage. Consider:

- protocol and expected RNA content;
- cell versus nucleus data;
- tissue-specific mitochondrial expression;
- whether a rare state is biologically plausible or damaged;
- whether retention differs systematically by condition or batch.

The barcode-rank figure is a diagnostic of the supplied cell matrix. It is not a replacement for upstream cell calling on an unfiltered droplet matrix.

### Doublet detection

```yaml
doublets:
  method: "scrublet"       # or "scdblfinder"
  expected_rate: 0.06
```

Run doublet detection per capture. Set the expected rate using the loading concentration and chemistry documentation, not a universal constant. The workflow passes this rate to Scrublet or to scDblFinder's `dbr` parameter. Compare predicted rates across samples and inspect whether flagged cells express incompatible lineage programs.

Neither method is ground truth. A large deviation can reflect incorrect expected rate, strong biological heterogeneity, too few cells, or upstream problems.

### Normalization

```yaml
normalization:
  method: "log"            # or "scran"
  target_sum: 10000
```

Log normalization is a practical default for visualization and marker ranking. Scran estimates pooling-based size factors and uses the R environment. Raw counts remain separate for count-based models.

Never use normalized or scaled values for pseudobulk count modeling.

### Highly variable genes and PCA

```yaml
hvg_pca:
  n_top_genes: 2000
  flavor: "seurat_v3"
  n_pcs: 50
  n_pcs_use: 30
```

`seurat_v3` expects raw counts for HVG selection. `n_pcs` is the number computed; `n_pcs_use` is the number used downstream.

Use the elbow plot and PCA loadings as diagnostics. There is no universal requirement to explain 80% of the variance. Check whether selected PCs are dominated by mitochondrial, ribosomal, cell-cycle, stress, sex-linked, or technical genes.

### Batch correction

```yaml
batch:
  batch_key: ""            # blank means no correction
  method: "harmony"        # harmony, scvi, or combat
  n_latent: 30
```

Leave `batch_key` blank unless you have a defensible technical covariate with multiple levels.

| Situation | Action |
|---|---|
| One technical batch | Do not correct |
| Multiple balanced technical batches with shared biology | Consider correction and compare before/after diagnostics |
| Condition perfectly confounded with batch | Do not claim correction can separate them; redesign or limit the question |
| A biological variable is being called “batch” | Stop and reconsider the analysis goal |

Harmony corrects a PCA representation. scVI learns a latent representation from counts. ComBat modifies expression under a linear model. These methods are not interchangeable, and stronger mixing is not automatically better.

### Neighbors and UMAP

```yaml
embedding:
  n_neighbors: 15
  metric: "euclidean"
  min_dist: 0.3
  spread: 1.0
  use_rep: ""
```

Leave `use_rep` blank for the representation selected by the batch configuration. Set it explicitly only when you know the relevant AnnData key and want strict validation.

UMAP is stochastic visualization of a neighbor graph. Do not interpret its axes, island area, or long-range distances as quantitative biology.

### Clustering: fixed biological choice, multiple diagnostics

```yaml
clustering:
  resolutions: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
  default_resolution: 0.6
```

The workflow calculates all configured candidate resolutions and reports cluster counts and silhouette values. It stores `default_resolution` as the primary `leiden` clustering.

Silhouette is a diagnostic, not an automatic cell-type selector. Choose a final resolution by jointly reviewing:

- marker coherence and known biology;
- stability across nearby resolutions;
- cluster size and whether tiny clusters are reproducible across samples;
- representation across samples, batches, and conditions;
- doublet, mitochondrial, stress, and cell-cycle signals;
- whether the scientific question needs broad cell classes or fine states.

If you revise the resolution, record the reason. Do not search many resolutions and report only the most convenient result.

### Cluster markers and pseudobulk DE

```yaml
markers:
  method: "wilcoxon"
  n_genes: 50
  logfc_min: 0.25
  pct_min: 0.1

  condition_key: "condition"
  group1: "control"
  group2: "treated"
  min_cells_per_sample: 10
  min_replicates_per_group: 2
  min_total_count: 10
  covariates: []
```

There are two different analyses here.

**Cluster markers** compare cells in one cluster with cells outside it. Use these rankings to characterize clusters. Because cells are not independent biological replicates and clustering used the same expression data, treat marker p-values as descriptive ranking aids.

**Condition DE** aggregates raw counts within each cluster and biological sample, then compares `group2` with `group1` using PyDESeq2. A positive log2 fold change means higher expression in `group2`, subject to the recorded contrast.

The workflow refuses or skips a cluster when guardrails fail, including insufficient cells per sample, insufficient replicates, a missing group, invalid counts, or a confounded design. Within an otherwise eligible cluster, `min_total_count` removes genes with too little summed pseudobulk signal to model. Review:

```text
results/merged/07_markers/pseudobulk_status.csv
```

On a completed marker rule, this table is part of the result. A skipped cluster is not a negative DE result. If every requested comparison is ineligible or fails, the rule stops and prints the status table to `logs/merged/07_markers.log`, because Snakemake removes declared outputs from failed jobs.

List sample-level covariates in `markers.covariates`. They must exist in the sample sheet, have valid variation, and leave a full-rank design. The workflow does not turn a confounded design into an estimable one.

### Annotation is opt-in

Annotation settings include:

```yaml
annotation:
  marker_set: ""           # "", "human_pbmc", or "mouse_brain"
  markers_json: ""         # optional custom {cell_type: [genes]} JSON
  min_marker_coverage: 0.5
  min_score_margin: 0.1
  celltypist: false
  celltypist_model: "Immune_All_Low.pkl"
  celltypist_min_confidence: 0.5
  manual_tsv: ""
```

The numerical values above illustrate the meaning of the controls; begin with the thresholds shipped in `config/config.yaml` and justify any change.

When no marker set, custom markers, CellTypist model, or manual mapping is selected, annotation remains `Unknown`. This is safer than assigning an unrelated reference label.

Configure one primary annotation route at a time unless you are deliberately comparing methods. If several are enabled, the final `cell_type` precedence is manual mapping, then CellTypist, then marker scoring; retain the method-specific columns for review.

Choose only one built-in marker set compatible with the species and tissue:

- `human_pbmc` is a limited immune reference, not a universal human atlas;
- `mouse_brain` is a limited brain reference, not a universal mouse atlas.

A custom marker JSON has this structure:

```json
{
  "Cell type A": ["GENE1", "GENE2", "GENE3"],
  "Cell type B": ["GENE4", "GENE5", "GENE6"]
}
```

`min_marker_coverage` determines whether each candidate label has enough of its markers in the dataset to be scored. Labels below this coverage are excluded; if no label is eligible, the workflow fails so that a species, gene-identifier, tissue, or marker-reference mismatch cannot pass silently. Among eligible labels, a cluster whose leading score is non-positive or whose lead over the runner-up is below `min_score_margin` remains `Unknown`.

CellTypist is also opt-in. Select a model whose species, tissue, assay, and label granularity match the data. The model may need network access on first use. The workflow constructs CellTypist's required log1p expression normalized to 10,000 counts per cell directly from the preserved raw-count layer, regardless of the main normalization method. Predictions below `celltypist_min_confidence` remain `Unknown`.

For manual annotation, provide a tab-separated cluster map:

```tsv
cluster	cell_type
0	Cell type A
1	Unknown
2	Cell type B
```

The manual file must contain exactly one row for every cluster. Use `Unknown` explicitly for unresolved clusters; a missing cluster or duplicate row is an input error. Manual annotation is not “ground truth” unless it is independently validated.

### Composition is descriptive

The workflow writes `composition_by_sample.csv` after annotation. It describes the fraction of retained cells assigned to each label in each sample. Those fractions can be affected by dissociation, capture, QC, doublet removal, and sampling depth.

Do not attach condition p-values to cell counts or percentages as if cells were independent replicates. Replicate-aware differential-abundance analysis is outside this workflow.

### PAGA and RNA velocity

```yaml
trajectory:
  run_paga: false
  loom_file: ""
  barcode_map_tsv: ""
  min_shared_cells: 100
```

Enable PAGA only when a continuous biological process is plausible and relevant. PAGA summarizes cluster connectivity; it does not establish direction, ancestry, or elapsed time.

scVelo additionally requires a loom file with valid `spliced` and `unspliced` count layers. The workflow enforces a barcode and gene contract:

- each velocity cell maps unambiguously to one merged analysis cell;
- loom and analysis barcodes are unique;
- genes are matched by unique identifiers and explicitly intersected/reordered;
- at least `trajectory.min_shared_cells` cells must remain after matching;
- sufficient shared genes and spliced/unspliced counts must exist.

When loom barcodes do not already equal the merged IDs, set `trajectory.barcode_map_tsv`. It must contain exactly these columns:

```tsv
analysis_barcode	loom_barcode
sample01_AAACCCAAGAAACACT-1	AAACCCAAGAAACACT:sample01
sample01_AAACCCAAGAAACCAT-1	AAACCCAAGAAACCAT:sample01
```

Both columns must be one-to-one and unique. `analysis_barcode` is the merged AnnData ID; `loom_barcode` is the original ID in the loom file.

RNA velocity is sensitive to protocol, kinetics, population selection, and model assumptions. Treat arrows and latent time as hypotheses to validate, not proof of lineage.

## 8. Plan and run the workflow

Run all commands from the repository root.

### Step 1: review inputs before computing

Confirm:

- every path exists and is readable;
- all matrices contain raw integer counts;
- all rows belong to one species and compatible reference;
- sample IDs are unique and safe;
- condition and batch are not confounded;
- FASTQ pairs are complete when present;
- annotation and trajectory are disabled unless their references are ready.

### Step 2: perform a dry-run

```bash
snakemake --use-conda --cores 8 --dry-run --printshellcmds
```

Read the output. A dry-run should show the intended samples, environments, and targets. Do not proceed past a missing-input or validation error.

### Step 3: optionally inspect the DAG

If Graphviz is installed:

```bash
snakemake --dag | dot -Tsvg > dag.svg
```

The raw-FASTQ branch and expression-analysis branch should appear as separate paths collected by the default targets. Raw-read QC must not feed or filter the count-matrix analysis.

### Step 4: run locally

```bash
snakemake --use-conda --cores 8 --rerun-incomplete --printshellcmds
```

Snakemake will create rule environments on their first use. The first run can therefore spend substantial time solving and downloading packages before analysis begins.

`--cores` limits concurrent scheduling; it does not guarantee that every library uses less memory. On a low-memory machine, reduce concurrency:

```bash
snakemake --use-conda --cores 2 --rerun-incomplete --printshellcmds
```

### Step 5: monitor without guessing

Use Snakemake's summary and the rule logs:

```bash
snakemake --summary
tail -n 50 logs/merged/06_cluster.log
tail -n 50 logs/merged/07_markers.log
```

The last command shown in the terminal is not always the failed command. Snakemake reports the rule and log path; inspect that log first.

### Step 6: resume

After fixing the cause of a failed job, rerun the same full command. Snakemake preserves complete outputs and resumes missing or incomplete work.

Use `snakemake --unlock` only after confirming that no Snakemake process is active and a stale lock is the actual problem.

### Step 7: confirm the declared run completed

```bash
snakemake --summary
ls -lh results/report.html results/run_config.yaml results/run_manifest.json
```

Confirm that the expected targets are complete and the three run-level files are nonempty. This verifies workflow completion, not biological validity; the next section explains the required review.

### Cluster execution

Snakemake 8 uses executor plugins. A typical SLURM pattern is:

```bash
snakemake --use-conda --cores 8 --jobs 100 \
  --executor cluster-generic \
  --cluster-generic-submit-cmd \
  'sbatch --cpus-per-task={threads} --mem={resources.mem_mb}M'
```

Site policies, storage, modules, and resource profiles differ. Test a dry-run and a small target with your cluster administrator before launching a full study.

Each analysis rule declares a conservative thread count and `mem_mb` request. You can override a default in `config/config.yaml`, for example:

```yaml
resources:
  qc: {threads: 2, mem_mb: 12000}
  markers: {threads: 8, mem_mb: 48000}
```

These values inform Snakemake scheduling and cluster submission. On a local machine they do not enforce a hard operating-system memory limit.

## 9. Interpret every analysis stage

Open `results/report.html`, but also inspect the stage-specific tables, figures, H5AD objects, logs, and manifest. The report is a navigation aid, not a substitute for review.

### 9.1 Input preflight

Confirm that the manifest lists the intended samples and inputs. Investigate every warning.

Stop if you find:

- normalized or fractional input values;
- mixed species or incompatible gene identifiers;
- duplicate samples or cell barcodes;
- a condition that occurs in only one technical batch;
- an unexpected input or reference path.

### 9.2 Per-sample cell QC

Key artifacts:

```text
results/<sample>/00_qc/qc_stats.json
results/<sample>/00_qc/figures/
results/<sample>/00_qc/qc_filtered.h5ad
```

Review counts, detected genes, mitochondrial fraction, and retention for every sample separately.

Ask:

- Were mitochondrial genes detected under the configured naming convention?
- Are thresholds supported by the observed distributions?
- Is a sample an outlier in depth or retention?
- Does filtering preferentially remove one condition or batch?
- Could a high-RNA population be incorrectly removed as a count outlier?
- Could a real low-RNA population be mistaken for damaged cells?

The knee plot describes barcodes already present in the input. It does not rerun empty-droplet detection.

### 9.3 Doublet detection

Key artifacts:

```text
results/<sample>/01_doublets/no_doublets.h5ad
results/<sample>/01_doublets/figures/01_doublet_score_hist.png
```

Compare the predicted doublet fraction with the expected loading rate and across captures. Review whether high-scoring cells co-express incompatible lineage programs or simply represent a legitimate transitional state.

An unexpectedly high rate is a warning to investigate; an apparently normal rate is not proof that all doublets were removed.

### 9.4 Merge and feature universe

Key artifact:

```text
results/merged/00_merged/merged.h5ad
```

Confirm the total cell count equals the retained per-sample counts and that genes use one consistent naming system. The merge requires identical gene-name sets across samples, permits different row order, and stops with missing/extra examples when the sets differ. Inspect sample representation before and after merging.

Structural zeros caused by incompatible references can look like biology. This workflow therefore does not perform an outer merge or impute absent features as zero. If samples were quantified against different annotations, resolve that upstream and regenerate comparable matrices.

### 9.5 Normalization

Key artifacts:

```text
results/merged/02_normalize/normalized.h5ad
results/merged/02_normalize/figures/
```

`X` contains normalized log expression for visualization and marker ranking. `layers["counts"]` retains raw counts for count-based models.

Review whether library-size distributions become comparable without erasing expected biological differences. Normalization does not remove batch effects and does not make cells independent replicates.

### 9.6 HVGs and PCA

Key artifacts:

```text
results/merged/03_hvg_pca/hvg_pca.h5ad
results/merged/03_hvg_pca/figures/01_hvg_dispersion.png
results/merged/03_hvg_pca/figures/02_pca_elbow.png
results/merged/03_hvg_pca/figures/03_pca_scatter_qc.png
results/merged/03_hvg_pca/figures/04_pca_loadings.png
```

Inspect which genes drive early PCs. A PC dominated by mitochondrial, ribosomal, stress, cell-cycle, sex-linked, or immunoglobulin genes may be biological, technical, or both. Decide based on the study, not on a generic blacklist.

Check that `n_pcs_use` captures stable structure without including a long tail of noise. Repeat the analysis with a nearby value when conclusions are sensitive.

### 9.7 Batch correction

Key artifacts:

```text
results/merged/04_batch/batch_corrected.h5ad
results/merged/04_batch/figures/
```

Compare before and after correction using several colorings:

- technical batch;
- sample;
- condition;
- broad cell identity;
- QC metrics.

Good correction reduces unwanted technical separation while preserving biological neighborhoods and markers. Complete batch mixing is not the goal when batches legitimately differ in cell composition.

Stop if correction mixes known distinct lineages, removes a condition-specific signal that the design can estimate, or creates structure driven by only one sample.

### 9.8 Neighbor graph and UMAP

Key artifacts:

```text
results/merged/05_embedding/embedded.h5ad
results/merged/05_embedding/figures/
```

Use UMAP to inspect local neighborhoods, sample coverage, and QC gradients. Do not interpret:

- axis values;
- island area as abundance;
- visual distance as a calibrated expression distance;
- a bridge as proof of transition; or
- a separated island as proof of a novel cell type.

Reproducible conclusions should also be visible in marker expression, the high-dimensional representation, and sample-level summaries.

### 9.9 Leiden clustering

Key artifacts:

```text
results/merged/06_cluster/clustered.h5ad
results/merged/06_cluster/clustering_summary.csv
results/merged/06_cluster/figures/
```

The `leiden` column corresponds to `clustering.default_resolution`. Candidate resolutions and silhouette values are sensitivity diagnostics.

For each reported cluster, ask:

- Does it have coherent positive and negative markers?
- Is it represented by more than one biological sample?
- Is it dominated by one technical batch?
- Is it driven by QC, doublet, stress, or cell-cycle signals?
- Does it remain recognizable at nearby resolutions?
- Is its size plausible, and is there enough replication for downstream tests?

A high silhouette score can prefer a coarse partition. It does not define the biologically correct number of cell types.

### 9.10 Cluster markers

Key artifacts:

```text
results/merged/07_markers/cluster_markers.csv
results/merged/07_markers/figures/
```

Prioritize:

- direction and magnitude of log fold change;
- fraction expressing inside and outside the cluster;
- agreement among several markers;
- known negative markers;
- consistency across samples.

Do not identify a cluster from one famous gene alone. Ambient RNA, stress, doublets, and low-level contamination can produce isolated marker expression.

Cell-level adjusted p-values in cluster-marker tables are ranking aids, not evidence for a condition effect. Use the pseudobulk branch for replicate-level condition comparisons.

### 9.11 Pseudobulk condition DE

Key artifacts:

```text
results/merged/07_markers/pseudobulk_status.csv
results/merged/07_markers/DE/pseudobulk_cluster_*.csv
```

Read `pseudobulk_status.csv` first. For every cluster it should state whether analysis ran or was skipped and why.

If the marker rule stopped because no requested pseudobulk comparison completed, inspect `logs/merged/07_markers.log` instead; it retains the same per-cluster status information on that failure path.

For a completed contrast, verify:

- the biological sample counts in each group;
- cells contributing to each sample-level pseudobulk;
- the exact `group2` versus `group1` contrast;
- covariates and design rank;
- effect sizes and uncertainty;
- adjusted p-values, not raw p-values alone;
- whether results are driven by one replicate;
- whether adequate counts exist for the cell population.

A cluster skipped for insufficient cells or replicates has no supported DE conclusion. It is not evidence of “no difference.”

The workflow's two-group model is not a substitute for paired, longitudinal, nested, repeated-measure, interaction, or random-effect analysis when the experiment requires one.

### 9.12 Annotation

Key artifacts:

```text
results/merged/08_annotate/annotated.h5ad
results/merged/08_annotate/annotation_summary.csv
results/merged/08_annotate/marker_coverage.csv
results/merged/08_annotate/composition_by_sample.csv
results/merged/08_annotate/figures/
```

Review marker coverage, score margin, classifier confidence, and the number of cells and samples behind every label. `Unknown` is a valid and often responsible result.

Use an iterative process:

1. Inspect clusters and marker rankings.
2. Choose a biologically compatible reference or custom marker set.
3. Run annotation.
4. Review positive and negative markers, reference confidence, and sample consistency.
5. Add a manual mapping only where evidence is sufficient.
6. Leave unresolved or mixed clusters as `Unknown`, or perform a justified lineage-restricted reanalysis.

Do not rename a cluster only to match an expected story. Automated annotation transfers labels from a reference; it does not validate the reference's relevance.

### 9.13 Composition

`composition_by_sample.csv` is the correct starting point because it preserves the sample as the unit of review.

Check whether apparent condition differences are consistent across independent samples and whether QC retention could explain them. Do not test cell-level contingency tables as if every cell were an independent replicate.

Formal differential abundance requires a replicate-aware compositional or neighborhood model and an experimental design that supports it. That inference is outside this workflow.

### 9.14 PAGA and scVelo

Key artifacts, when enabled:

```text
results/merged/09_trajectory/trajectory.h5ad
results/merged/09_trajectory/figures/
```

For PAGA, confirm that the analyzed population represents a plausible continuous process. Connectivity between clusters is not direction.

For scVelo, confirm the barcode mapping, shared cells, gene intersection, spliced/unspliced layer integrity, and lineage restriction. Review velocity confidence and model diagnostics. Do not use a global velocity graph across unrelated lineages merely because the command runs.

### 9.15 HTML report and run manifest

Key artifacts:

```text
results/report.html
results/run_config.yaml
results/run_manifest.json
```

The report gathers figures and summary tables for navigation. Review the source tables before making quantitative claims.

`run_config.yaml` is the effective configuration seen by Snakemake after command-line or profile overrides. Use this snapshot, rather than an edited source config, to reconstruct what the completed run requested.

The manifest records the effective configuration snapshot and hash, sample-sheet hash, input identity, declared-environment hashes, detected report-environment package versions, Git SHA, seed, and warnings. Preserve it with the results. Confirm that it describes the run you are interpreting, especially after changing configuration or rerunning selected stages.

## 10. Understand the output files

```text
results/
├── report.html
├── run_config.yaml
├── run_manifest.json
├── multiqc/
│   └── multiqc_report.html
├── <sample>/
│   ├── raw_fastq_qc/                  optional
│   │   ├── fastqc/
│   │   └── fastq_screen/
│   ├── 00_qc/
│   │   ├── qc_filtered.h5ad
│   │   ├── qc_stats.json
│   │   └── figures/
│   └── 01_doublets/
│       ├── no_doublets.h5ad
│       └── figures/
└── merged/
    ├── 00_merged/merged.h5ad
    ├── 02_normalize/normalized.h5ad
    ├── 03_hvg_pca/hvg_pca.h5ad
    ├── 04_batch/batch_corrected.h5ad
    ├── 05_embedding/embedded.h5ad
    ├── 06_cluster/
    │   ├── clustered.h5ad
    │   └── clustering_summary.csv
    ├── 07_markers/
    │   ├── with_markers.h5ad
    │   ├── cluster_markers.csv
    │   ├── pseudobulk_status.csv
    │   └── DE/pseudobulk_cluster_*.csv
    ├── 08_annotate/
    │   ├── annotated.h5ad
    │   ├── annotation_summary.csv
    │   ├── marker_coverage.csv
    │   └── composition_by_sample.csv
    └── 09_trajectory/                  optional
        ├── trajectory.h5ad
        └── figures/

logs/
├── <sample>/
└── merged/
```

### Which H5AD should you use?

- Use `qc_filtered.h5ad` to inspect one sample after cell QC.
- Use `no_doublets.h5ad` to inspect one sample after doublet removal.
- Use `normalized.h5ad` for the first merged normalized state.
- Use `clustered.h5ad` for clustering diagnostics before annotation.
- Use `annotated.h5ad` as the main downstream object after annotation review.
- Use `trajectory.h5ad` for PAGA output. When scVelo is enabled, this object contains only the validated cell-and-gene intersection and includes the velocity layers and results.

Keep raw input matrices and the manifest. Intermediate H5AD files are useful for debugging but can consume substantial storage.

## 11. Troubleshoot safely

### Snakemake reports a missing input file

Check the exact `path` in `config/samples.tsv`, your current working directory, file permissions, and spelling. Run from the repository root. Do not create an empty placeholder matrix.

### Preflight says counts are fractional or negative

The input is probably normalized, scaled, or otherwise transformed. Return to the upstream count-producing workflow and export the original molecule counts. Do not round normalized values.

### Preflight reports mixed species or incompatible genes

Split species into separate runs. For same-species annotation differences, reconcile gene identifiers and reference versions upstream and document the mapping.

### No mitochondrial genes are detected

Compare `var_names` with `mt_prefix` and `species`. Ensembl IDs require an explicit annotation strategy. Do not interpret zero mitochondrial percentage until the naming mismatch is resolved.

### Nearly all cells fail QC

Inspect the unfiltered distributions and `qc_stats.json`. Common causes include normalized input, wrong species prefixes, inappropriate hard cutoffs, nuclei analyzed with cell defaults, or a genuinely poor library. Change thresholds only with a biological and technical rationale.

### One sample loses far more cells than the others

Review raw-read QC, upstream cell calling, depth, mitochondrial fraction, and doublet rate. Check whether the sample is tied to one condition or batch; differential retention can bias composition and downstream comparisons.

### Doublet detection fails on a small sample

The method may not have enough cells or components. Inspect the log and sample size. Do not concatenate unrelated captures before doublet detection merely to make the method run.

### PCA requests too many components

Small datasets can have fewer usable cells or HVGs than configured PCs. Reduce `hvg_pca.n_pcs` and keep `n_pcs_use` below it.

### Batch correction fails or the requested representation is missing

Confirm that `batch.batch_key` exists, has multiple levels, and matches the selected method. An explicitly requested `embedding.use_rep` should exist; a missing corrected representation is an error, not a reason to silently reinterpret an uncorrected PCA.

### Batch correction erases a condition signal

Check for batch-condition confounding and whether a biological variable was used as the batch key. Compare uncorrected and corrected representations. If the design is confounded, restrict the claim rather than tuning integration until a preferred picture appears.

### A cluster has no markers

Inspect its size, neighboring resolutions, count depth, and whether it is defined by QC or batch. Verify that marker filters are not too strict. Do not lower every threshold until a desired marker appears.

### Marker annotation fails the coverage check

The selected marker dictionary is not compatible enough with the dataset to score safely. Confirm the global species, gene-symbol convention, tissue, and marker-reference choice. Prefer a compatible custom dictionary or corrected identifiers over lowering the coverage threshold until an unrelated reference runs.

### Many annotations are `Unknown`

This can be correct. Check gene-symbol compatibility, marker coverage, reference relevance, score margin, CellTypist confidence, and whether the cluster is mixed. Add a better reference or a carefully reviewed manual mapping; do not globally weaken uncertainty thresholds just to eliminate `Unknown`.

### CellTypist cannot download a model

The first use may require network access. Cache the intended model in advance or disable CellTypist and use a local marker/manual workflow. Record the exact model name and version.

### Pseudobulk DE is skipped

For a completed marker rule, open `pseudobulk_status.csv`. If the rule failed because every requested comparison was ineligible, open `logs/merged/07_markers.log`. Typical reasons include:

- fewer than `min_replicates_per_group` samples;
- too few cells from a sample in that cluster;
- one comparison group absent;
- a missing or constant covariate;
- a confounded or rank-deficient design;
- invalid raw counts.

Do not describe a skipped test as a nonsignificant result.

### A pseudobulk result is driven by one sample

Inspect sample-level counts and effect consistency. Consider whether the sample is a technical outlier, but do not remove it solely because it weakens significance. Use a prespecified, defensible exclusion rule and report sensitivity analyses.

### FastQ Screen cannot find a database

Copy the example configuration to the local configured path and replace every database placeholder with an existing Bowtie2 index prefix. Confirm `fastq_screen_reads` matches the assay read structure.

### MultiQC says no FASTQs were configured

That message is expected when both FASTQ columns are blank. It is a notice, not a passing raw-read QC result.

### scVelo finds too few shared cells

Check that `trajectory.barcode_map_tsv` has the exact `analysis_barcode` and `loom_barcode` columns, both are unique, and IDs refer to the intended samples. Confirm the merged sample prefix and any loom suffix conventions. Do not lower `min_shared_cells` until a mostly unmatched dataset runs.

### scVelo reports gene mismatch or missing layers

The loom file must contain unique gene identifiers plus `spliced` and `unspliced` count layers. Use the same reference annotation as the expression matrix. Do not assign layers by position without explicit gene matching.

### A job is killed or the machine runs out of memory

Reduce concurrent jobs with `--cores 1` or `--cores 2`. Inspect which rule failed. Large scran, ComBat, scVI, doublet, and H5AD operations may need a higher-memory machine. Reducing concurrency does not reduce the memory needed by one job.

### Conda cannot solve an environment

Retain the full solver output. Confirm channel order, platform support, free disk, and a current Conda/Mamba installation. Avoid repairing a rule environment with unrecorded `pip install` commands; update the environment definition and recreate it.

### Results look stale after a major config change

Use one clearly identified analysis directory per study version and retain its manifest. Snakemake tracks dependencies, but a major scientific change deserves a clean, separately named result set rather than silently overwriting the basis of an earlier interpretation.

## 12. Know the workflow's limits

The workflow does not provide:

- raw-read alignment or UMI quantification;
- cell calling from an unfiltered droplet matrix;
- ambient-RNA correction;
- chemistry-aware demultiplexing or donor assignment;
- sample-swap or genotype validation;
- guaranteed doublet truth;
- automatic selection of a biologically correct clustering resolution;
- universal cell-type annotation;
- replicate-aware differential-abundance inference;
- complex mixed-effect, longitudinal, nested, or interaction DE models;
- proof of lineage, direction, or developmental time;
- automatic biological validation;
- a substitute for an adequately replicated, balanced experimental design.

Other important considerations may include cell-cycle scoring, dissociation/stress programs, sex-linked expression, subclustering, reference mapping, pathway analysis, and independent validation. Add them only when they answer the study question and are documented prospectively.

The HTML report is not a publication-ready conclusion. It is a structured record for scientific review.

## 13. Final review checklist

Before accepting or sharing a result, confirm all of the following.

### Design

- [ ] `sample` represents the biological replicate for inference.
- [ ] Each condition has adequate independent replication.
- [ ] Conditions are balanced across technical batches.
- [ ] Covariates are justified and the design is not confounded.

### Inputs

- [ ] All matrices contain raw nonnegative integer counts in `X`.
- [ ] All samples use one species and compatible reference annotation.
- [ ] Gene and cell identifiers are unique and consistent.
- [ ] Upstream quantification, cell calling, and ambient-RNA decisions are recorded.

### QC and representation

- [ ] QC thresholds were reviewed per sample.
- [ ] Retention and doublet rates were compared across samples and conditions.
- [ ] HVGs and PCA loadings were inspected for nuisance drivers.
- [ ] Batch correction, if used, preserved expected biology.
- [ ] UMAP was treated as visualization rather than evidence by itself.

### Clustering and annotation

- [ ] The configured clustering resolution was justified biologically.
- [ ] Nearby resolutions and sample representation were reviewed.
- [ ] Labels are supported by several positive and negative markers.
- [ ] Weak or incompatible assignments remain `Unknown`.

### Inference

- [ ] Cluster markers were not presented as replicate-level condition tests.
- [ ] `pseudobulk_status.csv` was reviewed for every cluster.
- [ ] DE effect sizes, replicate consistency, contrast direction, and FDR were checked.
- [ ] Composition remained descriptive unless analyzed separately with a replicate-aware method.
- [ ] Trajectory or velocity claims were restricted to biologically plausible, validated populations.

### Provenance

- [ ] `results/run_manifest.json` matches the interpreted run.
- [ ] The effective configuration, sample sheet, input fingerprints, Conda-spec hashes, detected report-environment versions, Git SHA, seed, and warnings are retained.
- [ ] No raw data or generated results are mistaken for source code.

## 14. Glossary

**Adjusted p-value / FDR**

A multiple-testing-adjusted measure used to control the expected false-discovery proportion. It does not measure effect size or biological importance.

**Ambient RNA**

RNA molecules present in the droplet suspension that can be counted in cells where the transcript was not genuinely expressed. Correction is upstream of this workflow.

**AnnData / H5AD**

A cells-by-genes data object and its HDF5-backed file format. `X` is the main matrix; `obs`, `var`, `layers`, `obsm`, and `uns` store metadata and analysis state.

**Batch**

A technical grouping such as preparation day or sequencing run. A biological variable should not be relabeled as batch merely to remove it.

**Biological replicate**

An independently sampled experimental unit. This, not an individual cell, is the unit for condition-level inference.

**Cell barcode**

A sequence identifying the droplet or cell from which molecules originated.

**Cell calling**

The upstream decision that separates cell-containing barcodes from empty or background droplets.

**Condition**

The biological group being compared, such as control and treated.

**Confounding**

A design problem in which two effects always occur together and therefore cannot be estimated separately.

**Count matrix**

A cells-by-genes matrix of nonnegative integer molecule counts.

**Covariate**

A sample-level variable included in a statistical design to account for a known source of variation.

**Doublet**

A captured droplet or partition containing two or more cells, producing a mixed expression profile.

**Gene marker**

A gene enriched in a population and used, with other evidence, to characterize that population.

**HVG**

A highly variable gene selected to represent informative expression variation in dimensionality reduction.

**Integration / batch correction**

A method intended to reduce technical differences among datasets while preserving relevant biology.

**k-nearest-neighbor graph**

A graph linking each cell to nearby cells in a selected high-dimensional representation. It underlies UMAP and Leiden clustering.

**Leiden clustering**

A graph-community algorithm used to partition the cell-neighbor graph at a configured resolution.

**Log fold change**

The logarithm of an expression ratio. For log2 fold change, `1` means a two-fold increase and `-1` means a two-fold decrease.

**MEX**

The Matrix Exchange text format commonly used by 10x, usually represented by matrix, barcode, and feature files.

**PAGA**

Partition-based graph abstraction, which summarizes connectivity among clusters without establishing temporal direction.

**PCA**

Principal component analysis, a linear reduction that represents major axes of variation.

**Pseudobulk**

Raw counts summed within a cell population and biological sample so that the sample, not each cell, is the inferential unit.

**RNA velocity**

A model of transcriptional dynamics based on spliced and unspliced counts. Its direction and latent time depend on biological and model assumptions.

**Technical replicate**

A repeated measurement of the same biological material. It can assess technical variability but does not add independent biological replication.

**UMAP**

A nonlinear two-dimensional visualization of a neighbor graph. Its global geometry is not a calibrated biological distance.

**UMI**

A unique molecular identifier used to reduce amplification duplicates and estimate captured molecule counts.

For deeper methodological background, see [Single-cell best practices](https://www.sc-best-practices.org/), the [Scanpy documentation](https://scanpy.readthedocs.io/), and the [Snakemake documentation](https://snakemake.readthedocs.io/).
