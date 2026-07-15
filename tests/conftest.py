"""Shared test setup and safe loaders for numbered pipeline scripts."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIRECTORY = REPOSITORY_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIRECTORY))


@pytest.fixture
def load_script_module():
    """Load a numbered command-line script without executing its main block."""

    def load(filename: str):
        path = SCRIPTS_DIRECTORY / filename
        module_name = f"pipeline_test_{path.stem.replace('-', '_')}"
        spec = importlib.util.spec_from_file_location(module_name, path)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Could not load test module from {path}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    return load
