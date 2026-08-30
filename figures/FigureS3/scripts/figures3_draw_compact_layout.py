# Submission-aligned copy generated 2026-08-18
# Submission figure: FigureS3
# Role: Intermediate compact 24-panel Figure S3 layout

from __future__ import annotations

import json
import shutil
import warnings
import os
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc


warnings.filterwarnings("ignore", category=FutureWarning)

REPO_ROOT = Path(os.environ.get("SPATIALCCCBENCH_REPO_ROOT", Path(__file__).resolve().parents[3]))
PROJECT_DIR = Path(os.environ.get("SPATIALCCCBENCH_DATA_ROOT", REPO_ROOT))
OUT_DIR = Path(os.environ.get("SPATIALCCCBENCH_OUTPUT_DIR", Path(__file__).resolve().parents[1] / "results" / "generated"))
CACHE_DIR = OUT_DIR / "h5ad_readable_cache"
PANEL_DIR = OUT_DIR / "individual_panels"
COMPOSITE_DIR = OUT_DIR / "composites"
TABLE_DIR = OUT_DIR / "tables"


DATASETS = [
    {
        "panel": "a",
        "name": "DLPFC",
        "gene_adata": PROJECT_DIR / "adata_auto_corr_DLPFC_1.h5ad",
        "lr_adata": PROJECT_DIR / "adata_2_DLPFC.h5ad",
        "ligand": "SPP1",
        "receptor": "ITGB1",
        "lr_pair": "SPP1_ITGB1",
        "spot_size": 1.5,
        "start_index": 1,
    },
    {
        "panel": "b",
        "name": "MF",
        "gene_adata": PROJECT_DIR / "adata_MF_1.h5ad",
        "lr_adata": PROJECT_DIR / "adata_MF_2.h5ad",
        "ligand": "Pnoc",
        "receptor": "Oprl1",
        "lr_pair": "Pnoc_Oprl1",
        "spot_size": 0.3,
        "start_index": 13,
    },
]

# Top to bottom matches the requested compact FigureS3-like layout:
# random-like, homogeneous/gradient, heterogeneous/edge.
PER_GENE_MODES = [
    (
        "random_like_normal",
        "Random-like / normal",
        "colocalization",
        lambda moran, geary: (moran > 0) & (geary > 0.5) & (geary < 1.5),
    ),
    (
        "homogeneous_gradient",
        "Homogeneous / gradient",
        "colocalization",
        lambda moran, geary: (moran > 0) & (geary > 0) & (geary <= 0.5),
    ),
    (
        "heterogeneous_edge",
        "Heterogeneous / edge",
        "upheaval",
        lambda moran, geary: (moran > 0) & (geary >= 1.5),
    ),
]

LR_GEARY_MODES = [
    ("lr_geary_0p5_1p5", "Gearys=(0.5,1.5)", lambda geary: (geary > 0.5) & (geary < 1.5)),
    ("lr_geary_0_0p5", "Gearys=(0,0.5]", lambda geary: (geary > 0) & (geary <= 0.5)),
    ("lr_geary_1p5_inf", "Gearys=[1.5,+inf)", lambda geary: geary >= 1.5),
]

CONTINUOUS_CMAP = "RdYlGn_r"

PER_GENE_COLORS = {
    "interaction": "#1F77B4",
    "other": "#2CA02C",
    "ligand_upheaval": "#E6AB02",
    "receptor_upheaval": "#D62728",
    "ligand_colocalization": "#E6AB02",
    "receptor_colocalization": "#D62728",
}

LR_COLORS = {
    "Interaction": "#1F77B4",
    "Other": "#2CA02C",
    "Positive_Moran": "#E6AB02",
    "Geary": "#D62728",
}

RC_PARAMS = {
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "legend.title_fontsize": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.unicode_minus": False,
}


def ensure_dirs() -> None:
    for path in [CACHE_DIR, PANEL_DIR, COMPOSITE_DIR, TABLE_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def make_readable_h5ad(path: Path) -> Path:
    with h5py.File(path, "r") as h5:
        log1p_base = h5.get("uns/log1p/base")
        encoding_type = "" if log1p_base is None else log1p_base.attrs.get("encoding-type", "")
        if isinstance(encoding_type, bytes):
            encoding_type = encoding_type.decode("utf-8", errors="replace")
        has_legacy_null = str(encoding_type).lower() == "null"

    if has_legacy_null:
        fixed_path = CACHE_DIR / path.name
        if not fixed_path.exists() or fixed_path.stat().st_mtime < path.stat().st_mtime:
            shutil.copy2(path, fixed_path)
            with h5py.File(fixed_path, "r+") as h5:
                if "uns/log1p/base" in h5:
                    del h5["uns/log1p/base"]
        return fixed_path

    return path


def numeric_obs(adata: sc.AnnData, col: str) -> np.ndarray:
    return pd.to_numeric(adata.obs[col], errors="coerce").fillna(0).to_numpy()


def ensure_geary_log(adata: sc.AnnData, gene: str) -> str:
    log_col = f"{gene}_gearys(log)"
    if log_col in adata.obs.columns:
        return log_col
    geary_col = f"{gene}_gearys"
    adata.obs[log_col] = np.log10(np.clip(numeric_obs(adata, geary_col), 0, None) + 1)
    return log_col


def set_categorical(
    adata: sc.AnnData,
    col: str,
    values: np.ndarray,
    categories: list[str],
    palette: dict[str, str],
) -> None:
    adata.obs[col] = pd.Categorical(values, categories=categories)
    adata.uns[f"{col}_colors"] = [palette[cat] for cat in categories]


def add_right_legend(ax: plt.Axes, categories: list[str], palette: dict[str, str]) -> None:
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=palette[cat],
            markeredgecolor=palette[cat],
            markersize=3.0,
            label=cat,
        )
        for cat in categories
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.03, 0.5),
        frameon=False,
        handletextpad=0.25,
        borderaxespad=0,
        labelspacing=0.28,
        fontsize=6,
    )


def plot_continuous(
    adata: sc.AnnData,
    color: str,
    *,
    title: str,
    spot_size: float,
    ax: plt.Axes | None = None,
    output_path: Path | None = None,
    use_raw: bool | None = None,
) -> None:
    created = ax is None
    if created:
        fig, ax = plt.subplots(figsize=(2.65, 2.45))
    else:
        fig = ax.figure
    sc.pl.spatial(
        adata,
        color=color,
        title=title,
        size=spot_size,
        alpha_img=0.5,
        show=False,
        ax=ax,
        cmap=CONTINUOUS_CMAP,
        colorbar_loc="right",
        use_raw=use_raw,
    )
    if created and output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


def plot_categorical(
    adata: sc.AnnData,
    color: str,
    categories: list[str],
    palette: dict[str, str],
    *,
    title: str,
    spot_size: float,
    ax: plt.Axes | None = None,
    output_path: Path | None = None,
) -> None:
    created = ax is None
    if created:
        fig, ax = plt.subplots(figsize=(3.05, 2.45))
    else:
        fig = ax.figure
    sc.pl.spatial(
        adata,
        color=color,
        title=title,
        size=spot_size,
        alpha_img=0.5,
        show=False,
        ax=ax,
        legend_loc="none",
        palette=[palette[cat] for cat in categories],
    )
    add_right_legend(ax, categories, palette)
    if created and output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


def make_per_gene_interaction(
    adata: sc.AnnData,
    ligand: str,
    receptor: str,
    mode_key: str,
    prefix: str,
    mask_fn,
) -> tuple[str, list[str], dict[str, str], dict[str, int]]:
    ligand_mask = mask_fn(numeric_obs(adata, f"{ligand}_moran"), numeric_obs(adata, f"{ligand}_gearys"))
    receptor_mask = mask_fn(numeric_obs(adata, f"{receptor}_moran"), numeric_obs(adata, f"{receptor}_gearys"))
    interaction = ligand_mask & receptor_mask

    ligand_label = f"ligand_{prefix}"
    receptor_label = f"receptor_{prefix}"
    labels = np.full(adata.n_obs, "other", dtype=object)
    labels[ligand_mask] = ligand_label
    labels[receptor_mask] = receptor_label
    labels[interaction] = "interaction"

    categories = ["interaction", ligand_label, "other", receptor_label]
    palette = {cat: PER_GENE_COLORS[cat] for cat in categories}
    col = f"{ligand}_{receptor}_{mode_key}_per_gene_interaction"
    set_categorical(adata, col, labels, categories, palette)
    return col, categories, palette, {
        "ligand_count": int(ligand_mask.sum()),
        "receptor_count": int(receptor_mask.sum()),
        "interaction_count": int(interaction.sum()),
        "other_count": int((labels == "other").sum()),
    }


def make_lr_interaction(
    adata: sc.AnnData,
    lr_pair: str,
    mode_key: str,
    geary_label: str,
    mask_fn,
) -> tuple[str, list[str], dict[str, str], dict[str, int]]:
    moran = numeric_obs(adata, f"{lr_pair}_moran")
    geary = numeric_obs(adata, f"{lr_pair}_gearys")
    moran_mask = moran > 0
    geary_mask = mask_fn(geary)
    interaction = moran_mask & geary_mask

    labels = np.full(adata.n_obs, "Other", dtype=object)
    labels[moran_mask] = "Positive_Moran"
    labels[geary_mask] = geary_label
    labels[interaction] = "Interaction"

    categories = [geary_label, "Interaction", "Other", "Positive_Moran"]
    palette = {
        geary_label: LR_COLORS["Geary"],
        "Interaction": LR_COLORS["Interaction"],
        "Other": LR_COLORS["Other"],
        "Positive_Moran": LR_COLORS["Positive_Moran"],
    }
    col = f"{lr_pair}_{mode_key}_lr_integrated_interaction"
    set_categorical(adata, col, labels, categories, palette)
    return col, categories, palette, {
        "moran_positive_count": int(moran_mask.sum()),
        "geary_count": int(geary_mask.sum()),
        "interaction_count": int(interaction.sum()),
        "other_count": int((labels == "Other").sum()),
        "positive_moran_only_count": int((labels == "Positive_Moran").sum()),
        "geary_only_count": int((labels == geary_label).sum()),
    }


def check_columns(gene_adata: sc.AnnData, lr_adata: sc.AnnData, ligand: str, receptor: str, lr_pair: str) -> None:
    for gene in [ligand, receptor]:
        if gene not in gene_adata.var_names and (gene_adata.raw is None or gene not in gene_adata.raw.var_names):
            raise ValueError(f"Missing expression for {gene}")
        for suffix in ["moran", "gearys"]:
            if f"{gene}_{suffix}" not in gene_adata.obs.columns:
                raise ValueError(f"Missing {gene}_{suffix}")
    for suffix in ["moran", "gearys"]:
        if f"{lr_pair}_{suffix}" not in lr_adata.obs.columns:
            raise ValueError(f"Missing {lr_pair}_{suffix}")


def build_dataset_panels(config: dict[str, object]) -> tuple[list[list[dict[str, object]]], list[dict[str, object]]]:
    name = str(config["name"])
    ligand = str(config["ligand"])
    receptor = str(config["receptor"])
    lr_pair = str(config["lr_pair"])
    spot_size = float(config["spot_size"])

    gene_adata = sc.read_h5ad(make_readable_h5ad(Path(config["gene_adata"])), backed="r")
    lr_adata = sc.read_h5ad(make_readable_h5ad(Path(config["lr_adata"])), backed="r")
    check_columns(gene_adata, lr_adata, ligand, receptor, lr_pair)
    ligand_geary_log = ensure_geary_log(gene_adata, ligand)
    receptor_geary_log = ensure_geary_log(gene_adata, receptor)

    per_gene_specs = []
    count_rows = []
    for mode_key, mode_label, prefix, mask_fn in PER_GENE_MODES:
        col, categories, palette, counts = make_per_gene_interaction(
            gene_adata, ligand, receptor, mode_key, prefix, mask_fn
        )
        per_gene_specs.append(
            {
                "kind": "categorical",
                "adata": gene_adata,
                "color": col,
                "categories": categories,
                "palette": palette,
                "title": f"{lr_pair}_interaction",
                "spot_size": spot_size,
                "description": f"per_gene_{mode_key}",
            }
        )
        count_rows.append(
            {
                "dataset": name,
                "lr_pair": lr_pair,
                "panel_type": "per_gene_interaction",
                "mode": mode_label,
                **counts,
            }
        )

    lr_specs = []
    for mode_key, geary_label, mask_fn in LR_GEARY_MODES:
        col, categories, palette, counts = make_lr_interaction(lr_adata, lr_pair, mode_key, geary_label, mask_fn)
        lr_specs.append(
            {
                "kind": "categorical",
                "adata": lr_adata,
                "color": col,
                "categories": categories,
                "palette": palette,
                "title": f"{lr_pair}_interaction",
                "spot_size": spot_size,
                "description": f"lr_integrated_{mode_key}",
            }
        )
        count_rows.append(
            {
                "dataset": name,
                "lr_pair": lr_pair,
                "panel_type": "lr_integrated_interaction",
                "mode": geary_label,
                **counts,
            }
        )

    use_raw_ligand = True if gene_adata.raw is not None and ligand in gene_adata.raw.var_names else None
    use_raw_receptor = True if gene_adata.raw is not None and receptor in gene_adata.raw.var_names else None
    rows = [
        [
            {
                "kind": "continuous",
                "adata": gene_adata,
                "color": receptor,
                "title": receptor,
                "spot_size": spot_size,
                "use_raw": use_raw_receptor,
                "description": "receptor_raw_expression",
            },
            {
                "kind": "continuous",
                "adata": gene_adata,
                "color": ligand,
                "title": ligand,
                "spot_size": spot_size,
                "use_raw": use_raw_ligand,
                "description": "ligand_raw_expression",
            },
            per_gene_specs[0],
            lr_specs[0],
        ],
        [
            {
                "kind": "continuous",
                "adata": gene_adata,
                "color": ligand_geary_log,
                "title": ligand_geary_log,
                "spot_size": spot_size,
                "use_raw": None,
                "description": "ligand_local_geary_log",
            },
            {
                "kind": "continuous",
                "adata": gene_adata,
                "color": receptor_geary_log,
                "title": receptor_geary_log,
                "spot_size": spot_size,
                "use_raw": None,
                "description": "receptor_local_geary_log",
            },
            per_gene_specs[1],
            lr_specs[1],
        ],
        [
            {
                "kind": "continuous",
                "adata": gene_adata,
                "color": f"{ligand}_moran",
                "title": f"{ligand}_moran",
                "spot_size": spot_size,
                "use_raw": None,
                "description": "ligand_local_moran",
            },
            {
                "kind": "continuous",
                "adata": gene_adata,
                "color": f"{receptor}_moran",
                "title": f"{receptor}_moran",
                "spot_size": spot_size,
                "use_raw": None,
                "description": "receptor_local_moran",
            },
            per_gene_specs[2],
            lr_specs[2],
        ],
    ]
    return rows, count_rows


def render_spec(spec: dict[str, object], ax: plt.Axes | None = None, output_path: Path | None = None) -> None:
    if spec["kind"] == "continuous":
        plot_continuous(
            spec["adata"],
            str(spec["color"]),
            title=str(spec["title"]),
            spot_size=float(spec["spot_size"]),
            use_raw=spec["use_raw"],
            ax=ax,
            output_path=output_path,
        )
    else:
        plot_categorical(
            spec["adata"],
            str(spec["color"]),
            spec["categories"],
            spec["palette"],
            title=str(spec["title"]),
            spot_size=float(spec["spot_size"]),
            ax=ax,
            output_path=output_path,
        )


def save_individual_panels(config: dict[str, object], rows: list[list[dict[str, object]]]) -> list[dict[str, object]]:
    name = str(config["name"])
    lr_pair = str(config["lr_pair"])
    dataset_dir = PANEL_DIR / name
    dataset_dir.mkdir(parents=True, exist_ok=True)
    path_rows = []
    idx = int(config["start_index"])
    for row in rows:
        for spec in row:
            safe_desc = str(spec["description"]).replace("/", "_").replace(" ", "_")
            output_path = dataset_dir / f"{idx:02d}_{name}_{lr_pair}_{safe_desc}.svg"
            render_spec(spec, output_path=output_path)
            path_rows.append(
                {
                    "dataset": name,
                    "lr_pair": lr_pair,
                    "panel_index": idx,
                    "panel_type": spec["description"],
                    "output_path": output_path.relative_to(OUT_DIR).as_posix(),
                }
            )
            idx += 1
    return path_rows


def save_dataset_composite(config: dict[str, object], rows: list[list[dict[str, object]]]) -> None:
    name = str(config["name"])
    lr_pair = str(config["lr_pair"])
    fig = plt.figure(figsize=(12.0, 7.15))
    gs = GridSpec(3, 4, figure=fig, wspace=0.28, hspace=0.30)
    for r in range(3):
        for c in range(4):
            render_spec(rows[r][c], ax=fig.add_subplot(gs[r, c]))
    fig.subplots_adjust(left=0.045, right=0.965, top=0.965, bottom=0.045)
    fig.savefig(COMPOSITE_DIR / f"{name}_{lr_pair}_12_panels_compact.svg", dpi=300, bbox_inches="tight")
    fig.savefig(COMPOSITE_DIR / f"{name}_{lr_pair}_12_panels_compact.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_full_composite(dataset_rows: list[tuple[dict[str, object], list[list[dict[str, object]]]]]) -> None:
    fig = plt.figure(figsize=(12.2, 14.25))
    gs = GridSpec(
        7,
        4,
        figure=fig,
        height_ratios=[1, 1, 1, 0.14, 1, 1, 1],
        wspace=0.28,
        hspace=0.30,
    )
    row_offsets = [0, 4]
    for dataset_idx, (config, rows) in enumerate(dataset_rows):
        offset = row_offsets[dataset_idx]
        fig.text(0.014, 0.982 if dataset_idx == 0 else 0.498, str(config["panel"]), fontsize=22, fontfamily="Arial")
        for r in range(3):
            for c in range(4):
                render_spec(rows[r][c], ax=fig.add_subplot(gs[offset + r, c]))
    fig.subplots_adjust(left=0.045, right=0.965, top=0.982, bottom=0.025)
    fig.savefig(COMPOSITE_DIR / "FigureS3_like_24_compact_layout.svg", dpi=300, bbox_inches="tight")
    fig.savefig(COMPOSITE_DIR / "FigureS3_like_24_compact_layout.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ensure_dirs()
    plt.rcParams.update(RC_PARAMS)

    all_count_rows = []
    all_path_rows = []
    dataset_rows = []
    for config in DATASETS:
        rows, count_rows = build_dataset_panels(config)
        dataset_rows.append((config, rows))
        all_count_rows.extend(count_rows)
        all_path_rows.extend(save_individual_panels(config, rows))
        save_dataset_composite(config, rows)
    save_full_composite(dataset_rows)

    pd.DataFrame(all_count_rows).to_csv(TABLE_DIR / "interaction_counts.csv", index=False)
    pd.DataFrame(all_path_rows).to_csv(TABLE_DIR / "panel_paths.csv", index=False)
    (OUT_DIR / "run_config.json").write_text(
        json.dumps(
            {
                "layout": "Two stacked panels. Within each panel, columns 1-2 are expression/local Geary/local Moran and columns 3-4 are per-gene/LR-integrated interactions.",
                "row_order": ["random-like / normal", "homogeneous / gradient", "heterogeneous / edge"],
                "DLPFC": "SPP1-ITGB1 only",
                "MF": "Pnoc-Oprl1 only",
                "spot_size": {"DLPFC": 1.5, "MF": 0.3},
                "font": "Arial",
                "svg_fonttype": "none",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Saved compact FigureS3-like layout to {OUT_DIR}")


if __name__ == "__main__":
    main()
