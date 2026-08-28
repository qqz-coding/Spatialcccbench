# Submission-aligned copy generated 2026-08-18
# Submission figure: Figure2
# Role: Figure 2C-D and 2I-J Method 1 autocorrelation distributions and tool proportions

from __future__ import annotations

import ast
import json
import os
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator, interp1d
from scipy.stats import gaussian_kde


warnings.filterwarnings("ignore")

REPO_ROOT = Path(os.environ.get("SPATIALCCCBENCH_REPO_ROOT", Path(__file__).resolve().parents[3]))
PROJECT_DIR = Path(os.environ.get("SPATIALCCCBENCH_DATA_ROOT", REPO_ROOT))
OUT_DIR = Path(os.environ.get("SPATIALCCCBENCH_OUTPUT_DIR", Path(__file__).resolve().parents[1] / "results" / "generated"))
FIG_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"
EXTRACTED_DIR = OUT_DIR / "extracted_tool_results"

LIGAND_CSV = PROJECT_DIR / "result_l_1_full.csv"
RECEPTOR_CSV = PROJECT_DIR / "result_r_1_full.csv"
LR_REF_CSV = PROJECT_DIR / "dataset" / "cellchatDB_human_sl.csv"
SPOT_INFO_CSV = PROJECT_DIR / "dataset" / "lymph_spot_info.csv"

ANALYSIS_DATASET = "st_lymphnode"
GINI_THRESHOLD = 0.4

TOOL_LIST = [
    "Squidpy",
    "CellAgentChat",
    "SpaTalk",
    "SpatialDM",
    "stLearn",
    "stLearn_without_spotmixture",
    "COMMOT",
    "Giotto",
    "Baseline_2",
    "Baseline_1",
]

COUNT_COLUMNS = [
    "gradient_count",
    "edge_count",
    "normal_count",
    "weak_count",
]

COUNT_LABELS = {
    "gradient_count": "Homogeneous interaction",
    "edge_count": "Heterogenous interaction",
    "normal_count": "Random-like interaction",
    "weak_count": "Weak auto correlated interaction",
}

COUNT_COLORS = {
    "gradient_count": "pink",
    "edge_count": "red",
    "normal_count": "lightblue",
    "weak_count": "#FFEA00",
}

TOOL_BAR_COLUMNS = [
    "Homogeneous interaction",
    "Random-like interaction",
    "Heterogenous interaction",
    "Weak auto correlated interaction",
    "notfind",
]

TOOL_BAR_COLORS = ["pink", "lightblue", "red", "#FFEA00", "#90EE90"]


def ensure_dirs() -> None:
    for path in (FIG_DIR, TABLE_DIR, EXTRACTED_DIR):
        path.mkdir(parents=True, exist_ok=True)


def clean_gene_name(value: object) -> str:
    return str(value).strip().strip("'\"")


def parse_array(value: object) -> np.ndarray:
    if isinstance(value, np.ndarray):
        return value.astype(float)
    if isinstance(value, list):
        return np.asarray(value, dtype=float)
    text = str(value).strip()
    try:
        parsed = ast.literal_eval(text)
    except (SyntaxError, ValueError):
        parsed = ast.literal_eval(
            text.replace("nan", "None").replace("NaN", "None")
        )
    arr = np.asarray(parsed, dtype=float)
    return arr


def load_auto_corr(path: Path, gene_col: str) -> pd.DataFrame:
    usecols = [
        gene_col,
        "gini_index",
        "morani",
        "morani_p",
        "local_morani",
        "local_gearysC",
    ]
    df = pd.read_csv(path, usecols=usecols)
    df[gene_col] = df[gene_col].map(clean_gene_name)
    for col in ("local_morani", "local_gearysC"):
        df[col] = df[col].map(parse_array)
    return df


def validate_vector_lengths(result_ligand: pd.DataFrame, result_receptor: pd.DataFrame) -> int:
    lengths = set()
    for df in (result_ligand, result_receptor):
        for col in ("local_morani", "local_gearysC"):
            lengths.update(df[col].map(len).unique().tolist())
    if len(lengths) != 1:
        raise ValueError(f"Local arrays have inconsistent lengths: {sorted(lengths)}")
    return int(next(iter(lengths)))


def build_df_all(result_ligand: pd.DataFrame, result_receptor: pd.DataFrame) -> pd.DataFrame:
    ligand_names = set(result_ligand["ligand"])
    receptor_names = set(result_receptor["receptor"])

    lr_ref = pd.read_csv(LR_REF_CSV, index_col=0)
    lr_ref["ligand"] = lr_ref["ligand"].map(clean_gene_name)
    lr_ref["receptor"] = lr_ref["receptor"].map(clean_gene_name)
    lr_ref = lr_ref[
        lr_ref["ligand"].isin(ligand_names)
        & lr_ref["receptor"].isin(receptor_names)
    ][["ligand", "receptor"]].drop_duplicates()

    ligand_meta = result_ligand.set_index("ligand")
    receptor_meta = result_receptor.set_index("receptor")

    wide_ligands = set(
        result_ligand.loc[
            (result_ligand["gini_index"] > GINI_THRESHOLD)
            & (result_ligand["morani"] > 0),
            "ligand",
        ]
    )
    low_ligands = set(
        result_ligand.loc[
            (result_ligand["morani"] <= 0)
            | (result_ligand["gini_index"] <= GINI_THRESHOLD),
            "ligand",
        ]
    )
    wide_receptors = set(
        result_receptor.loc[
            (result_receptor["gini_index"] > GINI_THRESHOLD)
            & (result_receptor["morani"] > 0),
            "receptor",
        ]
    )
    low_receptors = set(
        result_receptor.loc[
            (result_receptor["morani"] <= 0)
            | (result_receptor["gini_index"] <= GINI_THRESHOLD),
            "receptor",
        ]
    )

    rows = []
    for ligand, receptor in lr_ref[["ligand", "receptor"]].itertuples(index=False):
        l_moran = ligand_meta.at[ligand, "local_morani"]
        l_geary = ligand_meta.at[ligand, "local_gearysC"]
        r_moran = receptor_meta.at[receptor, "local_morani"]
        r_geary = receptor_meta.at[receptor, "local_gearysC"]

        edge_count = gradient_count = normal_count = weak_count = 0

        if ligand in wide_ligands or receptor in wide_receptors:
            edge = ((l_moran > 0) & (l_geary >= 1.5)) | (
                (r_moran > 0) & (r_geary >= 1.5)
            )
            gradient = (
                (l_moran > 0)
                & (l_geary > 0)
                & (l_geary <= 0.5)
            ) | (
                (r_moran >= 0)
                & (r_geary > 0)
                & (r_geary <= 0.5)
            )
            normal = (
                (l_moran > 0)
                & (l_geary > 0.5)
                & (l_geary < 1.5)
            ) | (
                (r_moran > 0)
                & (r_geary > 0.5)
                & (r_geary < 1.5)
            )
            edge_count = int(edge.sum())
            gradient_count = int(gradient.sum())
            normal_count = int(normal.sum())

        if ligand in low_ligands or receptor in low_receptors:
            weak = ((l_moran <= 0) | (l_geary == 0)) & (
                (r_moran < 0) | (r_geary == 0)
            )
            weak_count = int(weak.sum())

        rows.append(
            {
                "LR": f"{ligand}_{receptor}",
                "gradient_count": gradient_count,
                "edge_count": edge_count,
                "normal_count": normal_count,
                "weak_count": weak_count,
            }
        )

    df_all = pd.DataFrame(rows).set_index("LR")[COUNT_COLUMNS].astype(float)
    zero_mask = (df_all[COUNT_COLUMNS] == 0).all(axis=1)
    df_all.loc[zero_mask, "weak_count"] = 1
    return df_all


def extract_tool_results() -> dict[str, pd.DataFrame]:
    sys.path.insert(0, str(REPO_ROOT))
    old_cwd = Path.cwd()
    os.chdir(PROJECT_DIR)
    try:
        from toolkit.preprocess import extract_result

        result_df_dict = extract_result(
            TOOL_LIST,
            ANALYSIS_DATASET,
            spot_info=str(SPOT_INFO_CSV),
        )
    finally:
        os.chdir(old_cwd)

    for tool, df in result_df_dict.items():
        df.to_csv(EXTRACTED_DIR / f"{tool}.csv", index=False)
    return result_df_dict


def build_grouped_df(df_all: pd.DataFrame) -> pd.DataFrame:
    df_reset = df_all.reset_index()
    grouped_frames = []
    group_specs = [
        ("gradient_count", "Local-co-localization"),
        ("edge_count", "Local-diffusion"),
        ("normal_count", "Global-co-localization"),
        ("weak_count", "Weak-spatial-auto-correlated"),
    ]
    for count_col, group_name in group_specs:
        group_df = df_reset.loc[df_reset[count_col] > 0].copy()
        group_df["group"] = group_name
        grouped_frames.append(group_df)
    return pd.concat(grouped_frames, ignore_index=True)


def calculate_tool_counts(
    df_all: pd.DataFrame, result_df_dict: dict[str, pd.DataFrame]
) -> tuple[pd.DataFrame, int]:
    df_grouped = build_grouped_df(df_all)
    total_count = len(df_grouped)
    rows = []

    for tool in TOOL_LIST:
        result_df = result_df_dict[tool]
        lr_list = set(
            result_df["LR_pairs"]
            .dropna()
            .astype(str)
            .str.strip()
            .str.replace("-", "_", regex=False)
        )
        filtered_df = df_grouped[df_grouped["LR"].isin(lr_list)]

        weak_counts = int(
            (filtered_df["group"] == "Weak-spatial-auto-correlated").sum()
        )
        gradient_counts = int(
            (filtered_df["group"] == "Local-co-localization").sum()
        )
        normal_counts = int(
            (filtered_df["group"] == "Global-co-localization").sum()
        )
        edge_counts = int((filtered_df["group"] == "Local-diffusion").sum())
        not_found = int(
            total_count
            - (gradient_counts + normal_counts + edge_counts + weak_counts)
        )

        rows.append(
            {
                "tool": tool,
                "Homogeneous interaction": gradient_counts,
                "Random-like interaction": normal_counts,
                "Heterogenous interaction": edge_counts,
                "Weak auto correlated interaction": weak_counts,
                "notfind": not_found,
            }
        )

    counts_df = pd.DataFrame(rows)
    counts_df["tool"] = counts_df["tool"].replace(
        "stLearn_without_spotmixture", "stLearn*"
    )
    counts_df = counts_df.set_index("tool")
    return counts_df[TOOL_BAR_COLUMNS], total_count


def savefig(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "pdf", "svg"):
        fig.savefig(
            FIG_DIR / f"{stem}.{ext}",
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="white",
        )
    plt.close(fig)


def soften_axes(ax, grid_axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#BDBDBD")
        ax.spines[side].set_linewidth(0.5)
    ax.tick_params(color="#BDBDBD", width=0.5)
    ax.grid(True, axis=grid_axis, alpha=0.2, linestyle="-", linewidth=0.35)
    ax.set_axisbelow(True)
    ax.set_facecolor("white")


def plot_tool_bar(df_percentage: pd.DataFrame, stem: str, legend: bool = True) -> None:
    fig, ax = plt.subplots(figsize=(6.2 if legend else 5, 5))
    plot_columns = [col for col in TOOL_BAR_COLUMNS if col in df_percentage.columns]
    plot_colors = [TOOL_BAR_COLORS[TOOL_BAR_COLUMNS.index(col)] for col in plot_columns]

    df_percentage[plot_columns].plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=plot_colors,
        linewidth=0.5,
        width=0.85,
    )

    ax.set_xlabel("Analysis tools")
    ax.set_ylabel("Percentage distribution")
    ax.set_ylim(0, 1)
    ax.tick_params(axis="x", rotation=65, labelsize=7, pad=1)
    ax.tick_params(axis="y", labelsize=7, pad=1)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.0%}"))
    soften_axes(ax, grid_axis="y")

    existing_legend = ax.get_legend()
    if legend:
        if existing_legend is not None:
            existing_legend.remove()
        ax.legend(
            plot_columns,
            title="Interaction Categories",
            frameon=False,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=7,
            title_fontsize=8,
        )
    elif existing_legend is not None:
        existing_legend.remove()

    fig.patch.set_facecolor("white")
    fig.tight_layout(pad=0.5)
    savefig(fig, stem)


def unique_empirical_cdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    sorted_values = np.sort(values)
    y = np.arange(1, len(sorted_values) + 1) / len(sorted_values)
    unique_values = []
    unique_y = []
    for idx, value in enumerate(sorted_values):
        if idx == 0 or value != sorted_values[idx - 1]:
            unique_values.append(value)
            unique_y.append(y[idx])
        else:
            unique_y[-1] = y[idx]
    return np.asarray(unique_values), np.asarray(unique_y)


def plot_count_cdf(df_all: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for col in COUNT_COLUMNS:
        unique_values, unique_y = unique_empirical_cdf(df_all[col].to_numpy())
        if len(unique_values) > 1:
            xs = np.linspace(unique_values.min(), unique_values.max(), 5000)
            if len(unique_values) >= 4:
                ys = PchipInterpolator(unique_values, unique_y)(xs)
            else:
                ys = interp1d(
                    unique_values,
                    unique_y,
                    kind="linear",
                    fill_value="extrapolate",
                )(xs)
            ys = np.maximum.accumulate(np.clip(ys, 0, 1))
            ax.plot(
                xs,
                ys,
                linewidth=1.2,
                label=COUNT_LABELS[col],
                color=COUNT_COLORS[col],
            )
    ax.set_xlabel("Spot count")
    ax.set_ylabel("Cumulative probability")
    ax.set_ylim(0, 1.05)
    soften_axes(ax, grid_axis="both")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout(pad=0.45)
    savefig(fig, "Method_1_CDF_counts")


def kde_smooth_cdf(values: np.ndarray, n_points: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(np.unique(values)) < 2:
        return np.asarray([0.0, 1.0]), np.asarray([0.0, 1.0])
    xs = np.linspace(0, 1, n_points)
    kde = gaussian_kde(values, bw_method="scott")
    pdf = kde(xs)
    cdf = np.cumsum(pdf)
    cdf = (cdf - cdf.min()) / (cdf.max() - cdf.min())
    cdf = np.maximum.accumulate(np.clip(cdf, 0, 1))
    return xs, cdf


def plot_ratio_cdf(df_all: pd.DataFrame) -> None:
    row_sums = df_all[COUNT_COLUMNS].sum(axis=1).replace(0, np.nan)
    ratio_df = df_all[COUNT_COLUMNS].div(row_sums, axis=0).fillna(0)
    ratio_df.to_csv(TABLE_DIR / "df_all_spot_count_ratios.csv")

    fig, ax = plt.subplots(figsize=(5, 3.5))
    for col in COUNT_COLUMNS:
        xs, ys = kde_smooth_cdf(ratio_df[col].to_numpy())
        ax.plot(
            xs,
            ys,
            linewidth=1.2,
            label=COUNT_LABELS[col],
            color=COUNT_COLORS[col],
        )
    ax.set_xlabel("Spot count proportion in interaction regions")
    ax.set_ylabel("Cumulative probability")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    soften_axes(ax, grid_axis="both")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout(pad=0.45)
    savefig(fig, "Method_1_CDF_ratios")


def write_category_summary(df_all: pd.DataFrame) -> None:
    totals = df_all[COUNT_COLUMNS].sum(axis=0)
    summary = pd.DataFrame(
        [
            {
                "category": COUNT_LABELS[col],
                "count_column": col,
                "total_count": totals[col],
                "proportion": totals[col] / totals.sum(),
                "lr_pairs": len(df_all),
            }
            for col in COUNT_COLUMNS
        ]
    )
    summary.to_csv(TABLE_DIR / "df_all_plot_summary.csv", index=False)


def main() -> None:
    ensure_dirs()
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "legend.title_fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    print("Loading ligand/receptor autocorrelation CSV files...")
    result_ligand = load_auto_corr(LIGAND_CSV, "ligand")
    result_receptor = load_auto_corr(RECEPTOR_CSV, "receptor")
    spot_count = validate_vector_lengths(result_ligand, result_receptor)
    print(
        f"Loaded {len(result_ligand)} ligands and {len(result_receptor)} receptors; "
        f"local vector length = {spot_count}."
    )

    print("Building df_all from local Moran/Geary arrays...")
    df_all = build_df_all(result_ligand, result_receptor)
    df_all.to_csv(TABLE_DIR / "df_all_from_result_l_1_full_result_r_1_full.csv")
    write_category_summary(df_all)
    plot_count_cdf(df_all)
    plot_ratio_cdf(df_all)

    print("Re-extracting tool results from result/st_lymphnode...")
    result_df_dict = extract_tool_results()
    counts_df, total_count = calculate_tool_counts(df_all, result_df_dict)
    counts_df.to_csv(TABLE_DIR / "tool_result_counts.csv")

    all_percentages = counts_df.div(total_count, axis=1).sort_values(
        "notfind", ascending=False
    )
    all_percentages.to_csv(TABLE_DIR / "tool_all_interaction_proportions.csv")
    plot_tool_bar(
        all_percentages,
        "All_interaction_types_of_CCC_result_1_full",
        legend=True,
    )

    value_columns = TOOL_BAR_COLUMNS[:-1]
    row_sums = counts_df[value_columns].sum(axis=1).replace(0, np.nan)
    spatial_percentages = counts_df[value_columns].div(row_sums, axis=0).fillna(0)
    spatial_percentages = spatial_percentages.sort_values(
        "Weak auto correlated interaction", ascending=True
    )
    spatial_percentages.to_csv(
        TABLE_DIR / "tool_spatial_correlated_proportions.csv"
    )
    plot_tool_bar(
        spatial_percentages,
        "Spatial_correlated_interaction_types_of_CCC_result_1_full",
        legend=True,
    )

    config = {
        "ligand_csv": str(LIGAND_CSV),
        "receptor_csv": str(RECEPTOR_CSV),
        "gini_threshold": GINI_THRESHOLD,
        "analysis_dataset": ANALYSIS_DATASET,
        "spot_count": spot_count,
        "lr_pairs_in_df_all": int(len(df_all)),
        "grouped_interaction_rows_for_tool_denominator": int(total_count),
        "weak_rule": "Figure1-A notebook: ((L_moran <= 0 or L_geary == 0) and (R_moran < 0 or R_geary == 0))",
        "outputs": {
            "figures": str(FIG_DIR),
            "tables": str(TABLE_DIR),
            "extracted_tool_results": str(EXTRACTED_DIR),
        },
    }
    (OUT_DIR / "run_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    print(f"Done. Outputs saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
