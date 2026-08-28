# Submission-aligned copy generated 2026-08-18
# Submission figure: Figure3
# Role: Figure 3A-C boundary adaptation rerun

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

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


def main():
    all_result = {}
    for level in LEVEL_LIST:
        print(f"Loading {level}...")
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
