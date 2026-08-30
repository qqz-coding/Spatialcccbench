"""Shared input, output, and reporting helpers for tool runners."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

import pandas as pd


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--analysis_dataset", required=True, help="Dataset label used under output_root.")
    parser.add_argument("--adata_path", required=True, help="Input AnnData h5ad file.")
    parser.add_argument("--LR_ref_path", required=True, help="CellChat-style ligand-receptor CSV.")
    parser.add_argument("--output_root", default="result", help="Root directory for benchmark outputs.")
    parser.add_argument("--cluster_key", default="cell_type", help="Cell-type annotation key.")
    parser.add_argument("--random_seed", type=int, default=0)


def output_directory(args: argparse.Namespace, tool_subdir: str) -> Path:
    path = Path(args.output_root).expanduser().resolve() / args.analysis_dataset / tool_subdir
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_lr_reference(path: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"ligand.symbol", "receptor.symbol"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"LR reference is missing columns: {', '.join(missing)}")
    return frame


def get_cell_type_matrix(adata: Any, cluster_key: str) -> pd.DataFrame:
    if cluster_key in adata.uns and isinstance(adata.uns[cluster_key], pd.DataFrame):
        matrix = adata.uns[cluster_key].reindex(adata.obs_names).fillna(0)
        return matrix.astype(float)
    if cluster_key in adata.obs:
        return pd.get_dummies(adata.obs[cluster_key], dtype=float)
    raise KeyError(
        f"Cell-type annotation '{cluster_key}' is absent from adata.obs and no composition matrix "
        f"exists in adata.uns['{cluster_key}']."
    )


def filter_housekeeping_genes(adata: Any) -> Any:
    names = adata.var_names.astype(str)
    keep = ~names.str.startswith(("MT-", "RPS", "RPL"))
    keep &= ~names.str.contains(r"^HB[^(P)]", regex=True)
    return adata[:, keep].copy()


def peak_memory_mb() -> float | None:
    try:
        import resource

        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return value / (1024 * 1024) if sys.platform == "darwin" else value / 1024
    except (ImportError, AttributeError):
        try:
            import psutil

            return psutil.Process().memory_info().rss / (1024 * 1024)
        except ImportError:
            return None


def report_resources(runtime_seconds: float) -> None:
    print(f"runtime_seconds: {runtime_seconds:.3f}")
    memory = peak_memory_mb()
    if memory is not None:
        print(f"peak_memory_mb: {memory:.2f}")
