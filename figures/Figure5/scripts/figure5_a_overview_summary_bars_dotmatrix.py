# Submission-aligned copy generated 2026-08-18
# Submission figure: Figure5
# Role: Figure 5A overview framework summary bars and dot matrix

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


DATA_ROOT = Path(os.environ.get("SPATIALCCCBENCH_DATA_ROOT", Path(__file__).resolve().parents[3]))
OUT_DIR = Path(os.environ.get("SPATIALCCCBENCH_OUTPUT_DIR", Path(__file__).resolve().parents[1] / "results" / "generated"))
SUMMARY_CSV = Path(os.environ.get("SPATIALCCCBENCH_SUMMARY_CSV", DATA_ROOT / "result_summary.csv"))

DISPLAY_TOOL_ORDER = [
    "Baseline_1",
    "Baseline_2",
    "CellAgentChat",
    "COMMOT",
    "Giotto",
    "SpaTalk",
    "SpatialDM",
    "Squidpy",
    "stLearn",
    "stLearn*",
]


RC_PARAMS = {
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
}


def min_max_normalize(values):
    values = pd.Series(values, dtype=float)
    min_val = values.min()
    max_val = values.max()
    if np.isclose(max_val, min_val):
        return pd.Series(1.0, index=values.index)
    return (values - min_val) / (max_val - min_val) + 1


def generate_gradient(start_hex, end_hex, steps=10):
    def hex_to_rgb(hex_color):
        hex_color = hex_color.lstrip("#")
        return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))

    start = np.array(hex_to_rgb(start_hex))
    end = np.array(hex_to_rgb(end_hex))
    delta = (end - start) / (steps - 1)
    gradient = [tuple((start + delta * i).astype(int)) for i in range(steps)]
    return [f"#{r:02x}{g:02x}{b:02x}" for r, g, b in gradient]


def create_color_map(colors):
    return LinearSegmentedColormap.from_list("custom_map", colors)


def get_color_ranks(values):
    values = pd.Series(values, dtype=float)
    ranks = values.rank(ascending=True)
    rank_range = ranks.max() - ranks.min()
    if np.isclose(rank_range, 0):
        return pd.Series(0.5, index=values.index)
    normalized_ranks = (ranks - ranks.min()) / rank_range
    return 1 - np.power(normalized_ranks, 0.7)


def save_figure(fig, stem):
    for suffix in ("png", "pdf", "svg"):
        is_svg = suffix == "svg"
        fig.savefig(
            OUT_DIR / f"{stem}.{suffix}",
            dpi=300,
            bbox_inches="tight",
            transparent=is_svg,
            facecolor="none" if is_svg else "white",
            edgecolor="none" if is_svg else "white",
        )
    plt.close(fig)


def prepare_result():
    result = pd.read_csv(SUMMARY_CSV, index_col=0)
    result = result.reindex([tool for tool in DISPLAY_TOOL_ORDER if tool in result.index])
    result["coverage"] = np.log10(result["Count"])
    result["tool"] = result.index
    result["Spatial_auto"] = result[["Spatial_auto_1", "Spatial_auto_2"]].mean(axis=1)
    result["Spatial_auto"] = min_max_normalize(100 - result["Spatial_auto"])
    result["Enrichment"] = min_max_normalize(result["Enrichment"])
    if "time" in result.columns:
        result["time"] = np.log10(result["time"])
    if "memory" in result.columns:
        result["memory"] = np.log10(result["memory"])
    return result


def build_gradients():
    gradients = [
        generate_gradient("#e74c3c", "#3498db", 10),
        generate_gradient("#2ecc71", "#e67e22", 10),
        generate_gradient("#e74c3c", "#9b59b6", 10),
        generate_gradient("#e67e22", "#f1c40f", 10),
        generate_gradient("#008000", "#2ecc71", 10),
        generate_gradient("#000080", "#9b59b6", 10),
        generate_gradient("#000000", "#ecf0f1", 10),
        generate_gradient("#8B0000", "#FFB6C1", 10),
        generate_gradient("#3498db", "#000000", 10),
        generate_gradient("#f1c40f", "#3498db", 10),
    ]
    color_dict = {
        "coverage": gradients[3],
        "Precision": gradients[0],
        "Recall": gradients[0],
        "F1 Score": gradients[0],
        "Enrichment": gradients[4],
        "Variance level": gradients[2],
        "Spatial_auto": gradients[5],
        "time": gradients[6],
        "memory": gradients[7],
    }
    return gradients, color_dict


def available_metrics(result, metrics):
    return [metric for metric in metrics if metric[0] in result.columns]


def plot_bar_summary(result, color_dict):
    metrics = [
        ("coverage", "Coverage", create_color_map(color_dict["coverage"])),
        ("Precision", "Precision", create_color_map(color_dict["Precision"])),
        ("Recall", "Recall", create_color_map(color_dict["Recall"])),
        ("F1 Score", "F1 score", create_color_map(color_dict["F1 Score"])),
        ("Enrichment", "Biological feature", create_color_map(color_dict["Enrichment"])),
        ("Variance level", "Spatial variance", create_color_map(color_dict["Variance level"])),
        ("Spatial_auto", "Spatial auto corr", create_color_map(color_dict["Spatial_auto"])),
        ("time", "Time cost", create_color_map(color_dict["time"])),
        ("memory", "Memory cost", create_color_map(color_dict["memory"])),
    ]
    metrics = available_metrics(result, metrics)
    methods = result.index

    fig, axes = plt.subplots(1, len(metrics), figsize=(1.05 * len(metrics), len(methods) * 0.36 + 1.4))
    axes = np.atleast_1d(axes)
    for i, (col, title, cmap) in enumerate(metrics):
        ax = axes[i]
        values = result[col]
        colors = [cmap(x) for x in get_color_ranks(values)]
        ax.barh(range(len(methods)), values, color=colors, edgecolor="none", linewidth=0)
        ax.set_title(title, fontsize=8, pad=6, y=-0.36, rotation=60, ha="center", va="bottom")
        ax.set_yticks(range(len(methods)))
        ax.set_yticklabels(methods if i == 0 else [], fontsize=7)
        ax.invert_yaxis()
        ax.tick_params(axis="x", length=0, labelbottom=False)
        ax.tick_params(axis="y", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(wspace=0.12, bottom=0.24)
    save_figure(fig, "bar_summary")


def plot_dot_matrix(result, metrics, stem, width_per_col=0.72):
    methods = result.index
    fig, axes = plt.subplots(1, len(metrics), figsize=(width_per_col * len(metrics), len(methods) * 0.42 + 1.3))
    axes = np.atleast_1d(axes)
    size_scale = np.array([45, 90, 135, 180, 225, 270])

    for i, (col, title, cmap) in enumerate(metrics):
        ax = axes[i]
        values = result[col]
        colors = [cmap(x) for x in get_color_ranks(values)]
        if np.isclose(values.max(), values.min()):
            sizes = np.repeat(size_scale.mean(), len(values))
        else:
            sizes = np.interp(values, [values.min(), values.max()], [size_scale[0], size_scale[-1]])
        ax.scatter(
            x=[1] * len(methods),
            y=range(len(methods)),
            s=sizes,
            c=colors,
            edgecolor="none",
            linewidth=0,
            alpha=0.9,
        )
        ax.set_title(title, fontsize=8, pad=6, y=-0.36, rotation=60, ha="center", va="bottom")
        ax.set_xlim(0.5, 1.5)
        ax.set_yticks(range(len(methods)))
        ax.set_yticklabels(methods if i == 0 else [], fontsize=7)
        ax.invert_yaxis()
        ax.xaxis.set_visible(False)
        ax.tick_params(axis="y", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if i == 0:
            ax.spines["left"].set_visible(True)
            ax.spines["left"].set_linewidth(0.6)
        ax.set_facecolor("white")

    fig.patch.set_facecolor("white")
    fig.subplots_adjust(wspace=0.02, bottom=0.24)
    save_figure(fig, stem)


def main():
    plt.rcParams.update(RC_PARAMS)
    result = prepare_result()
    gradients, color_dict = build_gradients()

    plot_bar_summary(result, color_dict)

    noise_metrics = [
        ("overlap_F1", "Spot overlap", create_color_map(gradients[8])),
        ("dropout_F1", "Drop out", create_color_map(gradients[8])),
        ("non_specific_F1", "Non-specific", create_color_map(gradients[8])),
        ("up_F1", "Random\ninterpolation", create_color_map(gradients[8])),
        ("offset_F1", "Off set", create_color_map(gradients[8])),
        ("lack_F1", "Spot lack", create_color_map(gradients[8])),
        ("dropdown_edge", "Edge signal loss", create_color_map(gradients[9])),
        ("equal_edge", "Edge equalized", create_color_map(gradients[9])),
        ("gradient_edge", "Edge gradient change", create_color_map(gradients[9])),
    ]
    plot_dot_matrix(result, available_metrics(result, noise_metrics), "bar_noise_sum", width_per_col=0.7)

    dot_metrics = [
        ("coverage", "Coverage", create_color_map(color_dict["coverage"])),
        ("Enrichment", "Biological feature", create_color_map(color_dict["Enrichment"])),
        ("Precision", "Precision", create_color_map(color_dict["Precision"])),
        ("Recall", "Recall", create_color_map(color_dict["Recall"])),
        ("F1 Score", "F1 score", create_color_map(color_dict["F1 Score"])),
        ("Variance level", "Spatial variance", create_color_map(color_dict["Variance level"])),
        ("Spatial_auto", "Spatial auto corr", create_color_map(color_dict["Spatial_auto"])),
    ]
    plot_dot_matrix(result, available_metrics(result, dot_metrics), "dot_sum", width_per_col=0.74)

    print(f"Saved summary figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
