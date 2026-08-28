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
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(
            OUT_DIR / f"{stem}.{suffix}",
            dpi=300,
            bbox_inches="tight",
            transparent=True,
        )
    plt.close(fig)


def prepare_result():
    result = pd.read_csv(SUMMARY_CSV, index_col=0)
    result["coverage"] = np.log10(result["Count"])
    result["tool"] = result.index
    result["Spatial_auto"] = result[["Spatial_auto_1", "Spatial_auto_2"]].mean(axis=1)
    result["Spatial_auto"] = min_max_normalize(100 - result["Spatial_auto"])
    result["Enrichment"] = min_max_normalize(result["Enrichment"])
    if "time" in result.columns:
        result["time"] = np.log10(result["time"])
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
        ("Variance level", "Spatial_variance", create_color_map(color_dict["Variance level"])),
        ("Spatial_auto", "Spatial_auto_corr", create_color_map(color_dict["Spatial_auto"])),
        ("time", "Time_cost", create_color_map(color_dict["time"])),
        ("memory", "Mermory_cost", create_color_map(color_dict["memory"])),
    ]
    metrics = available_metrics(result, metrics)
    methods = result.index

    fig, axes = plt.subplots(1, len(metrics), figsize=(30, len(methods) * 1.5), facecolor="none")
    fig.patch.set_alpha(0)
    axes = np.atleast_1d(axes)
    for ax in axes:
        ax.set_facecolor("none")
        for spine in ax.spines.values():
            spine.set_visible(False)
    for i, (col, title, cmap) in enumerate(metrics):
        ax = axes[i]
        values = result[col]
        colors = [cmap(x) for x in get_color_ranks(values)]
        ax.barh(range(len(methods)), values, color=colors, edgecolor="none", linewidth=0.5)
        ax.set_title(title, fontsize=12, pad=10, color="Black", loc="center", y=-0.3, rotation=60,
                     horizontalalignment="center", verticalalignment="bottom")
        ax.set_yticks(range(len(methods)))
        ax.set_yticklabels(methods if i == 0 else [], fontsize=10)
        ax.invert_yaxis()
    save_figure(fig, "bar_summary")


def plot_dot_matrix(result, metrics, stem, figsize):
    methods = result.index
    fig, axes = plt.subplots(1, len(metrics), figsize=figsize, facecolor="none")
    axes = np.atleast_1d(axes)
    for i, ax in enumerate(axes):
        ax.set_facecolor("none")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.xaxis.set_visible(False)
        if i == 0:
            ax.spines["left"].set_visible(True)
            ax.yaxis.set_visible(True)
        else:
            ax.yaxis.set_visible(False)

    size_scale = np.array([80, 160, 240, 320, 400, 480])

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
            linewidth=0.5,
            alpha=0.8,
        )
        ax.set_title(title, fontsize=12, pad=10, color="Black", loc="center", y=-0.3, rotation=60,
                     horizontalalignment="center", verticalalignment="bottom")
        ax.set_yticks(range(len(methods)))
        ax.set_yticklabels(methods if i == 0 else [], fontsize=10)
        ax.invert_yaxis()
    fig.subplots_adjust(wspace=0, hspace=0)
    save_figure(fig, stem)


def main():
    plt.rcdefaults()
    plt.rcParams["svg.fonttype"] = "path"
    result = prepare_result()
    gradients, color_dict = build_gradients()

    plot_bar_summary(result, color_dict)

    noise_metrics = [
        ("overlap_F1", "Spot overlap", create_color_map(gradients[8])),
        ("dropout_F1", "Drop out", create_color_map(gradients[8])),
        ("non_specific_F1", "Non-specific", create_color_map(gradients[8])),
        ("up_F1", "Random\ninterpolation", create_color_map(gradients[8])),
        ("offset_F1", "Off_set", create_color_map(gradients[8])),
        ("lack_F1", "Spot lack", create_color_map(gradients[8])),
        ("dropdown_edge", "Edge_signal_loss", create_color_map(gradients[9])),
        ("equal_edge", "Edge_equalized", create_color_map(gradients[9])),
        ("gradient_edge", "Edge_gradient_change", create_color_map(gradients[9])),
    ]
    plot_dot_matrix(result, available_metrics(result, noise_metrics), "bar_noise_sum", (10, len(result) * 1))

    dot_metrics = [
        ("coverage", "Coverage", create_color_map(color_dict["coverage"])),
        ("Enrichment", "Biological feature", create_color_map(color_dict["Enrichment"])),
        ("Precision", "Precision", create_color_map(color_dict["Precision"])),
        ("Recall", "Recall", create_color_map(color_dict["Recall"])),
        ("F1 Score", "F1 score", create_color_map(color_dict["F1 Score"])),
        ("Variance level", "Spatial_variance", create_color_map(color_dict["Variance level"])),
        ("Spatial_auto", "spatial_auto_corr", create_color_map(color_dict["Spatial_auto"])),
    ]
    plot_dot_matrix(result, available_metrics(result, dot_metrics), "dot_sum", (7, len(result) * 1))

    print(f"Saved summary figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
