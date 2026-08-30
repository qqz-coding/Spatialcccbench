from __future__ import annotations

import argparse
import os
import sys
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
from scipy.stats import ks_2samp


REPO_ROOT = Path(__file__).resolve().parents[3]
FIGURE3_SCRIPTS = REPO_ROOT / "figures" / "Figure3" / "scripts"
sys.path.insert(0, str(FIGURE3_SCRIPTS))

from figure3_preprocess_helper import extract_result  # noqa: E402


LEVELS = [
    "DLPFC",
    "DLPFC_equal",
    "DLPFC_REVERS",
    "DLPFC_down_02",
    "DLPFC_down_04",
    "DLPFC_down_06",
    "DLPFC_down_08",
]
DEFAULT_TOOLS = [
    "Squidpy",
    "Baseline_1",
    "Baseline_2",
    "COMMOT",
    "stLearn",
    "stLearn_without_spotmixture",
    "SpatialDM",
    "CellAgentChat",
    "SpaTalk",
    "Giotto",
]
DISPLAY_NAMES = {"stLearn_without_spotmixture": "stLearn*"}
DEFAULT_COUNTS = REPO_ROOT / "figures" / "Figure3" / "inputs" / "boundary_ligand_counts.csv"


@contextmanager
def working_directory(path: Path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def ligand_counts(frame: pd.DataFrame) -> pd.Series:
    pairs = frame["LR_pairs"].astype(str).str.replace("HLA_", "HLA-", regex=False)
    ligands = pairs.str.split("_", n=1, expand=True)[0]
    return ligands.value_counts().astype(float)


def signal_loss(original: pd.Series, dropout_counts: list[pd.Series]) -> float:
    differences = []
    for perturbed in dropout_counts:
        aligned = perturbed.reindex(original.index, fill_value=0)
        differences.append(float((original - aligned).mean()))
    return float(pd.Series(differences).mean())


def shared_ks(original: pd.Series, perturbed: pd.Series) -> float:
    shared = original.index.intersection(perturbed.index)
    original_shared = original.reindex(shared)
    perturbed_shared = perturbed.reindex(shared)
    return float(ks_2samp(original_shared, perturbed_shared).statistic)


def raw_ligand_counts(data_root: Path, tools: list[str]) -> dict[str, dict[str, pd.Series]]:
    with working_directory(data_root):
        extracted = {
            level: extract_result(
                tools,
                level,
                spot_info="./dataset/DLPFC_cell_info.csv",
            )
            for level in LEVELS
        }

    return {
        level: {tool: ligand_counts(extracted[level][tool]) for tool in tools}
        for level in LEVELS
    }


def compact_ligand_counts(path: Path, tools: list[str]) -> dict[str, dict[str, pd.Series]]:
    table = pd.read_csv(path)
    required = {"scenario", "tool", "ligand", "count"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    result: dict[str, dict[str, pd.Series]] = {}
    for level in LEVELS:
        result[level] = {}
        for tool in tools:
            subset = table[(table["scenario"] == level) & (table["tool"] == tool)]
            result[level][tool] = subset.set_index("ligand")["count"].astype(float)
    return result


def calculate_metrics(counts_by_level: dict[str, dict[str, pd.Series]], tools: list[str]) -> pd.DataFrame:
    rows = []
    for tool in tools:
        counts = {
            level: counts_by_level[level][tool]
            for level in LEVELS
        }
        rows.append(
            {
                "tool": DISPLAY_NAMES.get(tool, tool),
                "dropdown_edge": signal_loss(
                    counts["DLPFC"],
                    [counts[f"DLPFC_down_{level}"] for level in ("02", "04", "06", "08")],
                ),
                "equal_edge": shared_ks(counts["DLPFC"], counts["DLPFC_equal"]),
                "gradient_edge": shared_ks(counts["DLPFC"], counts["DLPFC_REVERS"]),
            }
        )
    return pd.DataFrame(rows).set_index("tool")


def update_summary(summary_path: Path, metrics: pd.DataFrame) -> None:
    summary = pd.read_csv(summary_path, index_col=0)
    missing = metrics.index.difference(summary.index)
    if not missing.empty:
        raise KeyError(f"Tools missing from summary: {missing.tolist()}")
    for column in ("dropdown_edge", "equal_edge", "gradient_edge"):
        summary.loc[metrics.index, column] = metrics[column]
    summary.to_csv(summary_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recompute Figure3/Figure5 boundary metrics from existing tool results."
    )
    parser.add_argument("--data-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--tools", nargs="+", default=DEFAULT_TOOLS)
    parser.add_argument("--counts-csv", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument(
        "--raw-results",
        action="store_true",
        help="Extract the full local result tree instead of using the compact latest-result table.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "figures" / "Figure5" / "inputs" / "boundary_metrics_recomputed.csv",
    )
    parser.add_argument("--update-summary", type=Path)
    args = parser.parse_args()

    if args.raw_results:
        counts_by_level = raw_ligand_counts(args.data_root.resolve(), args.tools)
    else:
        counts_by_level = compact_ligand_counts(args.counts_csv.resolve(), args.tools)
    metrics = calculate_metrics(counts_by_level, args.tools)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.output)
    if args.update_summary:
        update_summary(args.update_summary.resolve(), metrics)
    print(metrics.to_string())
    print(f"Saved: {args.output.resolve()}")


if __name__ == "__main__":
    main()
