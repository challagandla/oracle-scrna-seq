"""Cross-field preflight checks used while constructing the Snakemake DAG."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
REQUIRED_SAMPLE_COLUMNS = ("sample", "path", "species", "batch", "condition")


def _filled(value) -> bool:
    return pd.notna(value) and str(value).strip() not in {"", "nan", "None"}


def _require_file(path_value: str, label: str, errors: list[str]) -> None:
    if _filled(path_value) and not Path(str(path_value)).expanduser().exists():
        errors.append(f"{label} does not exist: {path_value}")


def _number(
    value,
    label: str,
    errors: list[str],
    *,
    minimum: float | None = None,
    maximum: float | None = None,
    integer: bool = False,
) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        errors.append(f"{label} must be numeric")
        return None
    if integer and not number.is_integer():
        errors.append(f"{label} must be an integer")
    if minimum is not None and number < minimum:
        errors.append(f"{label} must be at least {minimum:g}")
    if maximum is not None and number > maximum:
        errors.append(f"{label} must be at most {maximum:g}")
    return number


def _validate_count_input(path_value: str, label: str, errors: list[str]) -> None:
    if not _filled(path_value):
        return
    path = Path(str(path_value)).expanduser()
    if not path.exists():
        return
    if path.is_file() and path.suffix.lower() not in {".h5ad", ".h5"}:
        errors.append(f"{label} must be .h5ad, 10x .h5, or a 10x MEX directory")
    if path.is_dir():
        required_groups = (
            ("matrix.mtx", "matrix.mtx.gz"),
            ("barcodes.tsv", "barcodes.tsv.gz"),
            ("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"),
        )
        if any(not any((path / candidate).exists() for candidate in group) for group in required_groups):
            errors.append(f"{label} is not a complete 10x MEX directory: {path_value}")


def _perfectly_confounded(samples: pd.DataFrame, covariate: str, condition: str) -> bool:
    table = pd.crosstab(samples[covariate], samples[condition])
    return table.shape[1] > 1 and bool((table.gt(0).sum(axis=1) == 1).all())


def validate_project(config: dict, samples: pd.DataFrame) -> None:
    """Validate configuration, sample metadata, paths, and statistical design."""
    errors: list[str] = []

    if samples.empty:
        errors.append("samples_tsv must contain at least one sample row")

    missing = [column for column in REQUIRED_SAMPLE_COLUMNS if column not in samples.columns]
    if missing:
        errors.append(f"samples_tsv is missing required columns: {', '.join(missing)}")
    else:
        for column in REQUIRED_SAMPLE_COLUMNS:
            if samples[column].map(_filled).eq(False).any():
                errors.append(f"samples_tsv column '{column}' contains blank values")

    if "sample" in samples:
        identifiers = samples["sample"].astype(str)
        bad_ids = identifiers[~identifiers.map(lambda value: bool(SAMPLE_ID_RE.fullmatch(value)))]
        if not bad_ids.empty:
            errors.append(
                "sample IDs must match [A-Za-z0-9][A-Za-z0-9_.-]*; invalid: "
                + ", ".join(sorted(bad_ids.unique()))
            )
        duplicates = identifiers[identifiers.duplicated()].unique()
        if len(duplicates):
            errors.append("sample IDs must be unique; duplicates: " + ", ".join(duplicates))

    configured_species = str(config.get("species", "")).lower()
    if configured_species not in {"human", "mouse"}:
        errors.append("species must be 'human' or 'mouse'")
    if "species" in samples:
        observed_species = set(samples["species"].dropna().astype(str).str.lower())
        if observed_species - {"human", "mouse"}:
            errors.append("samples_tsv species values must be human or mouse")
        if observed_species and observed_species != {configured_species}:
            errors.append(
                f"all samples must match global species={configured_species}; found "
                + ", ".join(sorted(observed_species))
            )

    if "path" in samples:
        for row in samples.itertuples(index=False):
            label = f"count matrix for {row.sample}"
            _require_file(row.path, label, errors)
            _validate_count_input(row.path, label, errors)

    fastq_columns = {"raw_fastq_r1", "raw_fastq_r2"}
    present_fastq_columns = fastq_columns.intersection(samples.columns)
    if present_fastq_columns and present_fastq_columns != fastq_columns:
        errors.append("samples_tsv must contain both raw_fastq_r1 and raw_fastq_r2 columns")
    raw_samples: list[str] = []
    if fastq_columns.issubset(samples.columns):
        for row in samples.itertuples(index=False):
            r1, r2 = row.raw_fastq_r1, row.raw_fastq_r2
            if _filled(r1) != _filled(r2):
                errors.append(f"sample {row.sample} has only one raw FASTQ mate")
            if _filled(r1) and _filled(r2):
                raw_samples.append(str(row.sample))
                _require_file(r1, f"R1 FASTQ for {row.sample}", errors)
                _require_file(r2, f"R2 FASTQ for {row.sample}", errors)

    raw_cfg = config.get("raw_fastq_qc", {})
    screen_reads = raw_cfg.get("fastq_screen_reads", ["R2"])
    if not isinstance(screen_reads, list) or not screen_reads or not set(screen_reads) <= {"R1", "R2"}:
        errors.append("raw_fastq_qc.fastq_screen_reads must be a non-empty list of R1/R2")
    if raw_cfg.get("enabled", True) and raw_samples:
        screen_config = raw_cfg.get("fastq_screen_conf", "")
        _require_file(screen_config, "FastQ Screen configuration", errors)
        if str(screen_config).endswith(".example"):
            errors.append("copy the FastQ Screen example to a local configured file before raw QC")
        elif _filled(screen_config) and Path(str(screen_config)).exists():
            if "/path/to/" in Path(str(screen_config)).read_text(encoding="utf-8"):
                errors.append("FastQ Screen configuration still contains /path/to/ placeholders")

    _number(
        raw_cfg.get("fastq_screen_subset", 100000),
        "raw_fastq_qc.fastq_screen_subset",
        errors,
        minimum=1,
        integer=True,
    )

    _number(config.get("random_seed", 42), "random_seed", errors, minimum=0, integer=True)

    qc = config.get("qc", {})
    for key in ("min_cells", "min_counts", "min_genes"):
        _number(qc.get(key, 0), f"qc.{key}", errors, minimum=0, integer=True)
    for key in ("mad_counts", "mad_genes", "mad_mt"):
        _number(qc.get(key, 0), f"qc.{key}", errors, minimum=0)
    _number(qc.get("mt_hard", 0), "qc.mt_hard", errors, minimum=0, maximum=100)

    doublets = config.get("doublets", {})
    if doublets.get("method", "scrublet") not in {"scrublet", "scdblfinder"}:
        errors.append("doublets.method must be scrublet or scdblfinder")
    _number(
        doublets.get("expected_rate", 0.06),
        "doublets.expected_rate",
        errors,
        minimum=0,
        maximum=1,
    )

    normalization = config.get("normalization", {})
    if normalization.get("method", "log") not in {"log", "scran"}:
        errors.append("normalization.method must be log or scran")
    _number(normalization.get("target_sum", 10000), "normalization.target_sum", errors, minimum=1)

    hvg = config.get("hvg_pca", {})
    if hvg.get("flavor", "seurat_v3") not in {"seurat_v3", "seurat", "cell_ranger"}:
        errors.append("hvg_pca.flavor must be seurat_v3, seurat, or cell_ranger")
    _number(hvg.get("n_top_genes", 2000), "hvg_pca.n_top_genes", errors, minimum=2, integer=True)
    n_pcs = _number(hvg.get("n_pcs", 0), "hvg_pca.n_pcs", errors, minimum=2, integer=True)
    n_pcs_use = _number(
        hvg.get("n_pcs_use", 0), "hvg_pca.n_pcs_use", errors, minimum=2, integer=True
    )
    if n_pcs is not None and n_pcs_use is not None and n_pcs_use > n_pcs:
        errors.append("hvg_pca.n_pcs_use cannot exceed hvg_pca.n_pcs")

    batch = config.get("batch", {})
    batch_method = batch.get("method", "none")
    if batch_method not in {"none", "harmony", "scvi", "combat"}:
        errors.append("batch.method must be none, harmony, scvi, or combat")
    _number(batch.get("n_latent", 30), "batch.n_latent", errors, minimum=2, integer=True)
    batch_key = str(batch.get("batch_key") or "")
    if batch_key:
        if batch_key not in samples.columns:
            errors.append(f"batch.batch_key '{batch_key}' is not a samples_tsv column")
        elif samples[batch_key].nunique() < 2 and batch_method != "none":
            errors.append(f"batch correction needs at least two '{batch_key}' levels")
        elif (
            batch_method != "none"
            and "condition" in samples.columns
            and samples["condition"].nunique() > 1
            and _perfectly_confounded(samples, batch_key, "condition")
        ):
            errors.append(
                f"batch.batch_key '{batch_key}' is perfectly confounded with 'condition'; "
                "batch correction cannot recover this design"
            )

    embedding = config.get("embedding", {})
    representation = str(embedding.get("use_rep") or "")
    if representation not in {"", "X_pca", "X_pca_corrected", "X_scVI"}:
        errors.append("embedding.use_rep must be empty, X_pca, X_pca_corrected, or X_scVI")
    _number(embedding.get("n_neighbors", 15), "embedding.n_neighbors", errors, minimum=2, integer=True)
    min_dist = _number(embedding.get("min_dist", 0.3), "embedding.min_dist", errors, minimum=0)
    spread = _number(embedding.get("spread", 1.0), "embedding.spread", errors, minimum=0)
    if min_dist is not None and spread is not None and min_dist > spread:
        errors.append("embedding.min_dist cannot exceed embedding.spread")

    clustering = config.get("clustering", {})
    try:
        resolutions = [float(value) for value in clustering.get("resolutions", [])]
        default_resolution = float(clustering.get("default_resolution", -1))
    except (TypeError, ValueError):
        resolutions, default_resolution = [], -1
        errors.append("clustering resolutions must be numeric")
    if any(value <= 0 for value in resolutions) or len(set(resolutions)) != len(resolutions):
        errors.append("clustering.resolutions must contain unique positive values")
    if not resolutions or default_resolution not in resolutions:
        errors.append("clustering.default_resolution must appear in clustering.resolutions")

    markers = config.get("markers", {})
    if markers.get("method", "wilcoxon") not in {"wilcoxon", "t-test"}:
        errors.append("markers.method must be wilcoxon or t-test")
    _number(markers.get("n_genes", 50), "markers.n_genes", errors, minimum=1, integer=True)
    _number(markers.get("logfc_min", 0.25), "markers.logfc_min", errors, minimum=0)
    _number(markers.get("pct_min", 0.1), "markers.pct_min", errors, minimum=0, maximum=1)
    _number(
        markers.get("min_cells_per_sample", 10),
        "markers.min_cells_per_sample",
        errors,
        minimum=1,
        integer=True,
    )
    _number(
        markers.get("min_replicates_per_group", 2),
        "markers.min_replicates_per_group",
        errors,
        minimum=2,
        integer=True,
    )
    _number(
        markers.get("min_total_count", 10),
        "markers.min_total_count",
        errors,
        minimum=1,
        integer=True,
    )
    condition_key = str(markers.get("condition_key") or "")
    if condition_key:
        group1, group2 = str(markers.get("group1") or ""), str(markers.get("group2") or "")
        if condition_key not in samples.columns:
            errors.append(f"markers.condition_key '{condition_key}' is not a samples_tsv column")
        elif not group1 or not group2 or group1 == group2:
            errors.append("pseudobulk DE needs two distinct non-empty group labels")
        else:
            levels = set(samples[condition_key].astype(str))
            if not {group1, group2} <= levels:
                errors.append(f"pseudobulk groups {group1}/{group2} are not both present")
            minimum = int(markers.get("min_replicates_per_group", 2))
            counts = samples[condition_key].astype(str).value_counts()
            for group in (group1, group2):
                if int(counts.get(group, 0)) < minimum:
                    errors.append(f"pseudobulk group '{group}' has fewer than {minimum} samples")
        covariates = markers.get("covariates", []) or []
        if not isinstance(covariates, list):
            errors.append("markers.covariates must be a list of samples_tsv columns")
        else:
            for covariate in covariates:
                if covariate not in samples.columns:
                    errors.append(f"pseudobulk covariate '{covariate}' is not a samples_tsv column")
                elif condition_key in samples.columns and _perfectly_confounded(samples, covariate, condition_key):
                    errors.append(
                        f"pseudobulk covariate '{covariate}' is perfectly confounded with "
                        f"'{condition_key}'"
                    )

    annotation = config.get("annotation", {})
    if annotation.get("marker_set", "") not in {"", "human_pbmc", "mouse_brain"}:
        errors.append("annotation.marker_set must be empty, human_pbmc, or mouse_brain")
    marker_set = annotation.get("marker_set", "")
    expected_marker_set = {"human": "human_pbmc", "mouse": "mouse_brain"}.get(
        configured_species
    )
    if marker_set and expected_marker_set and marker_set != expected_marker_set:
        errors.append(
            f"annotation.marker_set '{marker_set}' is incompatible with "
            f"species='{configured_species}'"
        )
    if marker_set and _filled(annotation.get("markers_json", "")):
        errors.append("choose annotation.marker_set or annotation.markers_json, not both")
    _number(
        annotation.get("min_marker_coverage", 0.5),
        "annotation.min_marker_coverage",
        errors,
        minimum=0,
        maximum=1,
    )
    _number(
        annotation.get("min_score_margin", 0.1),
        "annotation.min_score_margin",
        errors,
        minimum=0,
    )
    _number(
        annotation.get("celltypist_min_confidence", 0.5),
        "annotation.celltypist_min_confidence",
        errors,
        minimum=0,
        maximum=1,
    )
    _require_file(annotation.get("markers_json", ""), "annotation markers JSON", errors)
    _require_file(annotation.get("manual_tsv", ""), "manual annotation TSV", errors)
    markers_json = annotation.get("markers_json", "")
    if _filled(markers_json) and Path(str(markers_json)).is_file():
        try:
            marker_dictionary = json.loads(Path(str(markers_json)).read_text(encoding="utf-8"))
            if not isinstance(marker_dictionary, dict) or not marker_dictionary:
                raise ValueError("it must be a non-empty object")
            clean_labels = set()
            for label, genes in marker_dictionary.items():
                clean_label = str(label).strip()
                clean_genes = {
                    str(gene).strip() for gene in genes
                } if isinstance(genes, list) else set()
                clean_genes.discard("")
                if not clean_label or len(clean_genes) < 2:
                    raise ValueError("every label needs at least two unique genes")
                if clean_label in clean_labels:
                    raise ValueError("marker labels must be unique after trimming")
                clean_labels.add(clean_label)
        except (OSError, json.JSONDecodeError, ValueError) as error:
            errors.append(f"annotation markers JSON is invalid: {error}")

    manual_tsv = annotation.get("manual_tsv", "")
    if _filled(manual_tsv) and Path(str(manual_tsv)).is_file():
        try:
            manual = pd.read_csv(manual_tsv, sep="\t", dtype=str)
            if set(manual.columns) != {"cluster", "cell_type"}:
                errors.append("manual annotation TSV must contain exactly cluster and cell_type")
            elif manual.empty:
                errors.append("manual annotation TSV must contain at least one mapping")
        except (OSError, pd.errors.ParserError) as error:
            errors.append(f"manual annotation TSV is invalid: {error}")

    trajectory = config.get("trajectory", {})
    loom_file = trajectory.get("loom_file", "")
    barcode_map = trajectory.get("barcode_map_tsv", "")
    _require_file(loom_file, "trajectory loom file", errors)
    _require_file(barcode_map, "trajectory barcode map", errors)
    if _filled(barcode_map):
        try:
            mapping = pd.read_csv(barcode_map, sep="\t", dtype=str)
            required_mapping = {"analysis_barcode", "loom_barcode"}
            if set(mapping.columns) != required_mapping:
                errors.append(
                    "trajectory barcode map must contain exactly analysis_barcode and "
                    "loom_barcode columns"
                )
            else:
                for column in required_mapping:
                    mapping[column] = mapping[column].str.strip()
                if (
                    mapping.empty
                    or mapping[list(required_mapping)].isna().any().any()
                    or mapping[list(required_mapping)].eq("").any().any()
                    or any(mapping[column].duplicated().any() for column in required_mapping)
                ):
                    errors.append("trajectory barcode map must be complete and one-to-one")
        except (OSError, pd.errors.ParserError) as error:
            errors.append(f"trajectory barcode map is invalid: {error}")
    if _filled(barcode_map) and not _filled(loom_file):
        errors.append("trajectory.barcode_map_tsv requires trajectory.loom_file")
    _number(
        trajectory.get("min_shared_cells", 100),
        "trajectory.min_shared_cells",
        errors,
        minimum=1,
        integer=True,
    )

    for rule_name, values in config.get("resources", {}).items():
        if rule_name not in {
            "raw_fastqc", "fastq_screen", "multiqc", "qc", "doublets", "merge",
            "normalize", "hvg_pca", "batch", "embedding", "cluster", "markers",
            "annotate", "trajectory", "report",
        }:
            errors.append(f"resources contains unknown rule '{rule_name}'")
            continue
        if not isinstance(values, dict):
            errors.append(f"resources.{rule_name} must be a mapping")
            continue
        for key in ("threads", "mem_mb"):
            if key in values:
                _number(
                    values[key],
                    f"resources.{rule_name}.{key}",
                    errors,
                    minimum=1,
                    integer=True,
                )

    if errors:
        details = "\n  - ".join(errors)
        raise ValueError(f"Project preflight failed:\n  - {details}")
