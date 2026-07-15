#!/usr/bin/env python3
"""Generate a self-contained report and machine-readable run manifest."""

from __future__ import annotations

import argparse
import base64
import hashlib
import importlib.metadata
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml


WARNINGS = [
    "This workflow starts from a cell-by-gene raw-count matrix. Raw-read QC does not perform alignment, UMI quantification, cell calling, or ambient-RNA correction.",
    "Cluster-marker p-values treat cells as observations and are descriptive aids for annotation, not condition-level biological inference.",
    "Automated cell-type labels are provisional. Unknown labels and disagreements require marker review and biological validation.",
    "Cell-type proportions are descriptive unless a separate replicate-aware differential-abundance model is used.",
    "PAGA and RNA velocity do not by themselves prove lineage direction or temporal causality.",
]

TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>scRNA-seq workflow report</title>
<style>
body {{font-family:Arial,sans-serif;margin:0;background:#f5f7fa;color:#1f2933}}
header {{background:#17324d;color:white;padding:22px 40px}} main {{max-width:1200px;margin:auto;padding:28px 20px}}
section {{background:white;border-radius:8px;padding:22px;margin-bottom:22px;box-shadow:0 1px 4px #0002}}
h2 {{border-bottom:2px solid #2f80a8;padding-bottom:6px}} table {{border-collapse:collapse;width:100%;font-size:.86rem}}
th {{background:#2f80a8;color:white;text-align:left}} th,td {{padding:7px;border-bottom:1px solid #e5e7eb}}
.grid {{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:16px}} .card img {{width:100%}}
.stats {{display:flex;gap:12px;flex-wrap:wrap}} .stat {{background:#edf6fa;border-left:4px solid #2f80a8;padding:10px 16px}}
.warning {{background:#fff5db;border-left:4px solid #d99000;padding:9px 12px;margin:8px 0}}
</style></head><body>
<header><h1>scRNA-seq matrix-to-report workflow</h1><p>Generated {date} · commit {commit}</p></header><main>
<section><h2>Run summary</h2><div class="stats">
<div class="stat"><strong>{n_cells}</strong><br>cells</div><div class="stat"><strong>{n_genes}</strong><br>genes</div>
<div class="stat"><strong>{n_samples}</strong><br>samples</div><div class="stat"><strong>{n_clusters}</strong><br>clusters</div>
</div></section>
<section><h2>Interpretation boundaries</h2>{warnings}</section>
<section><h2>Per-sample QC retention</h2>{qc_table}</section>
<section><h2>Sample design</h2>{sample_table}</section>
<section><h2>Clustering diagnostics</h2>{cluster_table}</section>
<section><h2>Annotation summary</h2>{annotation_table}</section>
<section><h2>Annotation marker coverage</h2>{coverage_table}</section>
<section><h2>Pseudobulk status</h2>{de_table}</section>
<section><h2>Top descriptive cluster markers</h2>{marker_table}</section>
<section><h2>Figures</h2><div class="grid">{figures}</div></section>
</main></body></html>"""


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path_value, hash_contents=False):
    path = Path(path_value).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Cannot fingerprint missing input: {path}")
    if path.is_file():
        stat = path.stat()
        return {
            "path": str(path),
            "type": "file",
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": sha256_file(path) if hash_contents else None,
        }
    files = sorted(item for item in path.rglob("*") if item.is_file())
    digest = hashlib.sha256() if hash_contents else None
    size = 0
    latest_mtime = 0
    for item in files:
        stat = item.stat()
        size += stat.st_size
        latest_mtime = max(latest_mtime, stat.st_mtime_ns)
        if digest:
            digest.update(str(item.relative_to(path)).encode())
            digest.update(sha256_file(item).encode())
    return {
        "path": str(path),
        "type": "directory",
        "files": len(files),
        "size": size,
        "mtime_ns": latest_mtime,
        "sha256": digest.hexdigest() if digest else None,
    }


def git_commit():
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def software_versions():
    packages = [
        "anndata", "scanpy", "numpy", "pandas", "scipy", "snakemake",
        "scrublet", "harmonypy", "scvi-tools", "pydeseq2", "celltypist", "scvelo",
    ]
    versions = {"python": sys.version.split()[0]}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def build_manifest(args, samples, adata):
    config = yaml.safe_load(Path(args.config_yaml).read_text(encoding="utf-8"))
    input_columns = [column for column in ("path", "raw_fastq_r1", "raw_fastq_r2") if column in samples]
    inputs = []
    for column in input_columns:
        for value in pd.unique(samples[column].fillna("")):
            if str(value).strip():
                record = fingerprint(value, hash_contents=args.hash_inputs)
                record["role"] = f"sample_sheet.{column}"
                inputs.append(record)

    configured_inputs = {
        "raw_fastq_qc.fastq_screen_conf": config.get("raw_fastq_qc", {}).get(
            "fastq_screen_conf", ""
        ),
        "annotation.markers_json": config.get("annotation", {}).get("markers_json", ""),
        "annotation.manual_tsv": config.get("annotation", {}).get("manual_tsv", ""),
        "trajectory.loom_file": config.get("trajectory", {}).get("loom_file", ""),
        "trajectory.barcode_map_tsv": config.get("trajectory", {}).get(
            "barcode_map_tsv", ""
        ),
    }
    for role, value in configured_inputs.items():
        if not str(value).strip():
            continue
        path = Path(value).expanduser()
        if path.exists():
            record = fingerprint(value, hash_contents=args.hash_inputs)
            record["role"] = role
        else:
            record = {
                "path": str(path.resolve()),
                "role": role,
                "exists": False,
            }
        inputs.append(record)
    return {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(),
        "config_sha256": sha256_file(args.config_yaml),
        "samples_sha256": sha256_file(args.samples_tsv),
        "config": config,
        "samples": samples.fillna("").to_dict(orient="records"),
        "inputs": inputs,
        "input_content_hashing": bool(args.hash_inputs),
        "report_environment_software": software_versions(),
        "environment_specs": {
            Path(path).name: fingerprint(path, hash_contents=True)
            for path in args.environment_specs
        },
        "summary": {
            "cells": int(adata.n_obs),
            "genes": int(adata.n_vars),
            "samples": int(samples["sample"].nunique()),
            "clusters": int(adata.obs["leiden"].nunique()) if "leiden" in adata.obs else None,
        },
        "warnings": WARNINGS,
    }


def dataframe_html(table, max_rows=50):
    if table is None or table.empty:
        return "<p>Not available or not requested.</p>"
    return table.head(max_rows).to_html(index=False, border=0, escape=True, float_format="{:.4f}".format)


def image_card(path, caption):
    with open(path, "rb") as handle:
        encoded = base64.b64encode(handle.read()).decode()
    extension = Path(path).suffix.lstrip(".")
    return f'<div class="card"><img src="data:image/{extension};base64,{encoded}" alt="figure"><p>{caption}</p></div>'


def collect_qc(samples, results_dir):
    rows = []
    for sample in samples["sample"].astype(str):
        path = Path(results_dir) / sample / "00_qc" / "qc_stats.json"
        if not path.exists():
            continue
        stats = json.loads(path.read_text(encoding="utf-8"))
        before, after = stats["before"], stats["after"]
        rows.append(
            {
                "sample": sample,
                "cells_before": before["n_cells"],
                "cells_after": after["n_cells"],
                "retained_pct": 100 * after["n_cells"] / max(before["n_cells"], 1),
                "median_counts_after": after["median_counts"],
                "median_genes_after": after["median_genes"],
                "mean_mt_pct_after": after["mean_pct_mt"],
            }
        )
    return pd.DataFrame(rows)


def main(args):
    annotated = ad.read_h5ad(args.annotated)
    samples = pd.read_csv(args.samples_tsv, sep="\t", dtype=str)
    clusters = pd.read_csv(args.clust_csv)
    markers = pd.read_csv(args.markers)
    annotation = pd.read_csv(args.annotation_summary)
    marker_coverage = pd.read_csv(args.marker_coverage)
    de_status = pd.read_csv(args.de_status)

    manifest = build_manifest(args, samples, annotated)
    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    Path(args.manifest).write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")

    top_markers = (
        markers.sort_values(["cluster", "padj"])
        .groupby("cluster", observed=True)
        .head(3)
        [["cluster", "gene", "log2FC", "padj"]]
        if not markers.empty
        else markers
    )
    figures = []
    for path in sorted(Path(args.results_dir).rglob("*.png")):
        figures.append(image_card(path, str(path.relative_to(args.results_dir))))

    commit = manifest["git_commit"]
    html = TEMPLATE.format(
        date=manifest["generated_utc"],
        commit=commit[:12] if commit != "unavailable" else commit,
        n_cells=f"{annotated.n_obs:,}",
        n_genes=f"{annotated.n_vars:,}",
        n_samples=samples["sample"].nunique(),
        n_clusters=annotated.obs["leiden"].nunique() if "leiden" in annotated.obs else "—",
        warnings="".join(f'<div class="warning">{warning}</div>' for warning in WARNINGS),
        qc_table=dataframe_html(collect_qc(samples, args.results_dir)),
        sample_table=dataframe_html(samples.drop(columns=["path", "raw_fastq_r1", "raw_fastq_r2"], errors="ignore")),
        cluster_table=dataframe_html(clusters),
        annotation_table=dataframe_html(annotation),
        coverage_table=dataframe_html(marker_coverage),
        de_table=dataframe_html(de_status),
        marker_table=dataframe_html(top_markers),
        figures="\n".join(figures),
    )
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(html, encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotated", required=True)
    parser.add_argument("--markers", required=True)
    parser.add_argument("--clust-csv", required=True)
    parser.add_argument("--annotation-summary", required=True)
    parser.add_argument("--marker-coverage", required=True)
    parser.add_argument("--de-status", required=True)
    parser.add_argument("--samples-tsv", required=True)
    parser.add_argument("--config-yaml", required=True)
    parser.add_argument("--environment-specs", nargs="+", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--hash-inputs", action="store_true")
    main(parser.parse_args())
