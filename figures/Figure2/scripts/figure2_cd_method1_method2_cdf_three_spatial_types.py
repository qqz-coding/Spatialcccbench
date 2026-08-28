# Submission-aligned copy generated 2026-08-18
# Submission figure: Figure2
# Role: Figure 2C-D CDF curves for three spatial autocorrelation classes

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


REPO_ROOT = Path(os.environ.get("SPATIALCCCBENCH_REPO_ROOT", Path(__file__).resolve().parents[3]))
PROJECT_DIR = Path(os.environ.get("SPATIALCCCBENCH_DATA_ROOT", REPO_ROOT))
OUT_DIR = Path(os.environ.get("SPATIALCCCBENCH_OUTPUT_DIR", Path(__file__).resolve().parents[1] / "results" / "generated"))
FIG_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"

INPUTS = {
    "Method_1": PROJECT_DIR
    / "figure1a_reproduce_full"
    / "tables"
    / "df_all_from_result_l_1_full_result_r_1_full.csv",
    "Method_2": PROJECT_DIR
    / "figure1b_reproduce_full"
    / "tables"
    / "df_all_from_result_auto_corr_2.csv",
}

COUNT_COLUMNS = ["gradient_count", "edge_count", "normal_count"]
COUNT_LABELS = {
    "gradient_count": "Homogeneous interaction",
    "edge_count": "Heterogenous interaction",
    "normal_count": "Random-like interaction",
}
COUNT_COLORS = {
    "gradient_count": "pink",
    "edge_count": "red",
    "normal_count": "lightblue",
}


def ensure_dirs() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)


def kde_smooth_count_cdf(values: np.ndarray, n_points: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    unique_values = np.unique(values)
    if len(unique_values) < 2:
        x0 = float(unique_values[0]) if len(unique_values) else 0.0
        return np.asarray([x0, x0 + 1.0]), np.asarray([0.0, 1.0])

    xs = np.linspace(values.min(), values.max(), n_points)
    kde = gaussian_kde(values, bw_method="scott")
    pdf = kde(xs)
    cdf = np.cumsum(pdf)
    cdf = (cdf - cdf.min()) / (cdf.max() - cdf.min())
    cdf = np.maximum.accumulate(np.clip(cdf, 0, 1))
    return xs, cdf


def soften_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.grid(True, axis="both", color="#E0E0E0", linewidth=0.5, alpha=0.7)


def savefig(fig: plt.Figure, stem: str) -> None:
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(FIG_DIR / f"{stem}.{ext}", dpi=300, bbox_inches="tight")


def plot_cdf(df_all: pd.DataFrame, method_label: str) -> pd.DataFrame:
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["font.size"] = 10
    plt.rcParams["axes.labelsize"] = 11
    plt.rcParams["axes.titlesize"] = 12
    plt.rcParams["legend.fontsize"] = 9

    fig, ax = plt.subplots(figsize=(10, 6))
    rows = []
    for col in COUNT_COLUMNS:
        values = df_all[col].to_numpy(dtype=float)
        xs, ys = kde_smooth_count_cdf(values)
        ax.plot(
            xs,
            ys,
            linewidth=3,
            label=COUNT_LABELS[col],
            color=COUNT_COLORS[col],
            alpha=0.9,
        )
        rows.append(
            {
                "method": method_label,
                "count_column": col,
                "label": COUNT_LABELS[col],
                "color": COUNT_COLORS[col],
                "lr_pairs": int(len(values)),
                "nonzero_lr_pairs": int((values > 0).sum()),
                "zero_lr_pairs": int((values == 0).sum()),
                "min": float(np.nanmin(values)),
                "max": float(np.nanmax(values)),
                "total_spot_count": float(np.nansum(values)),
                "included_in_cdf": True,
            }
        )

    ax.set_xlabel("Value", fontsize=12, fontweight="bold")
    ax.set_ylabel("Cumulative Probability", fontsize=12, fontweight="bold")
    ax.set_title("KDE-Smoothed CDF Comparison", fontsize=14, fontweight="bold", pad=20)
    ax.legend(fontsize=11, frameon=False, loc="lower right")
    ax.set_ylim(0, 1.05)
    soften_axes(ax)
    fig.tight_layout()
    savefig(fig, f"{method_label}_CDF_three_spatial_autocorr_types")
    plt.close(fig)
    return pd.DataFrame(rows)


def run() -> None:
    ensure_dirs()
    all_rows = []
    for method_label, csv_path in INPUTS.items():
        df_all = pd.read_csv(csv_path)
        missing = [col for col in COUNT_COLUMNS if col not in df_all.columns]
        if missing:
            raise ValueError(f"Missing columns in {csv_path}: {missing}")
        all_rows.append(plot_cdf(df_all, method_label))

    summary = pd.concat(all_rows, ignore_index=True)
    summary.to_csv(TABLE_DIR / "three_spatial_autocorr_cdf_summary.csv", index=False)
    (OUT_DIR / "run_config.json").write_text(
        json.dumps(
            {
                "inputs": {method: str(path) for method, path in INPUTS.items()},
                "outputs": str(FIG_DIR),
                "cdf_method": "KDE-smoothed count CDF",
                "included_columns": COUNT_COLUMNS,
                "excluded_columns": ["weak_count"],
                "reason": "Only the three spatially autocorrelated interaction classes are drawn.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    run()
