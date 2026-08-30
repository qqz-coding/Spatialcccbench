# Submission-aligned copy generated 2026-08-18
# Submission figure: Figure3
# Role: Figure 3A-C boundary adaptation rerun

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

from figure3_boundary_adaptation_helpers import (
    analysis_edge_adoption_dropout,
    analysis_edge_adoption_equal,
    analysis_edge_adoption_gradient_change,
    deep_analysis,
)
from figure3_preprocess_helper import extract_result


TOOL_LIST = [
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

LEVEL_LIST = [
    "DLPFC",
    "DLPFC_equal",
    "DLPFC_REVERS",
    "DLPFC_down_02",
    "DLPFC_down_04",
    "DLPFC_down_06",
    "DLPFC_down_08",
]

DEFAULT_COUNTS = Path(__file__).resolve().parents[1] / "inputs" / "boundary_ligand_counts.csv"


def load_compact_counts(path: Path):
    counts = pd.read_csv(path)
    required = {"scenario", "tool", "ligand", "count"}
    missing = required.difference(counts.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")

    all_result = {}
    for level in LEVEL_LIST:
        all_result[level] = {}
        for tool in TOOL_LIST:
            subset = counts[(counts["scenario"] == level) & (counts["tool"] == tool)]
            repeated = subset.loc[subset.index.repeat(subset["count"].astype(int))]
            all_result[level][tool] = pd.DataFrame(
                {
                    "LR_pairs": repeated["ligand"].astype(str).map(lambda ligand: f"{ligand}_R").to_numpy(),
                    "cell_pairs": "",
                }
            )
    return all_result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts-csv", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--raw-results", action="store_true")
    args = parser.parse_args()

    if args.counts_csv.exists() and not args.raw_results:
        print(f"Loading compact latest-result counts: {args.counts_csv}")
        all_result = load_compact_counts(args.counts_csv)
    else:
        all_result = {}
        for level in LEVEL_LIST:
            print(f"Loading raw results for {level}...")
            all_result[level] = extract_result(
                TOOL_LIST,
                level,
                spot_info="./dataset/DLPFC_cell_info.csv",
            )

    for tool in TOOL_LIST:
        print(f"Plotting dropout: {tool}")
        analysis_edge_adoption_dropout(tool, all_result)
        plt.close("all")

    for tool in TOOL_LIST:
        print(f"Plotting edge equalization: {tool}")
        analysis_edge_adoption_equal(tool, all_result)
        plt.close("all")

    for tool in TOOL_LIST:
        print(f"Plotting gradient change: {tool}")
        analysis_edge_adoption_gradient_change(tool, all_result)
        plt.close("all")

    for analysis_target in ["DLPFC_equal", "DLPFC_REVERS"]:
        for tool in TOOL_LIST:
            print(f"Plotting deep analysis: {tool} / {analysis_target}")
            deep_analysis(all_result, tool, analysis_target=analysis_target)
            plt.close("all")


if __name__ == "__main__":
    main()
