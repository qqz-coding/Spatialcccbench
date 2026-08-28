# Submission-aligned copy generated 2026-08-18
# Submission figure: Figure2
# Role: Figure 2A-B CXCL10-CXCR3 local Moran and Geary spatial maps

from __future__ import annotations

import ast
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


REPO_ROOT = Path(os.environ.get("SPATIALCCCBENCH_REPO_ROOT", Path(__file__).resolve().parents[3]))
PROJECT_DIR = Path(os.environ.get("SPATIALCCCBENCH_DATA_ROOT", REPO_ROOT))
OUT_DIR = Path(os.environ.get("SPATIALCCCBENCH_OUTPUT_DIR", Path(__file__).resolve().parents[1] / "results" / "generated"))
SVG_DIR = OUT_DIR / "svg"
TABLE_DIR = OUT_DIR / "tables"

ADATA_PATH = PROJECT_DIR / "adata_1.h5ad"
LR_AUTO_CSV = PROJECT_DIR / "result_auto_corr_2.csv"

LIGAND = "CXCL10"
RECEPTOR = "CXCR3"
DISPLAY_PAIR = "CXCR3_CXCL10"
SOURCE_PAIR = f"{LIGAND}_{RECEPTOR}"

RC_PARAMS = {
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "legend.title_fontsize": 8,
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
    "axes.unicode_minus": False,
}

STATE_SPECS = [
    (
        "homogeneous_gradient",
        "Homogeneous / gradient",
        lambda moran, geary: (moran > 0) & (geary > 0) & (geary <= 0.5),
    ),
    (
        "random_like_normal",
        "Random-like / normal",
        lambda moran, geary: (moran > 0) & (geary > 0.5) & (geary < 1.5),
    ),
    (
        "heterogeneous_edge",
        "Heterogeneous / edge",
        lambda moran, geary: (moran > 0) & (geary > 1.5),
    ),
]


def ensure_dirs() -> None:
    SVG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)


def parse_numeric_list(value: object) -> np.ndarray:
    if isinstance(value, np.ndarray):
        return value.astype(float)
    if isinstance(value, list):
        return np.asarray(value, dtype=float)

    text = str(value).strip()
    try:
        return np.asarray(ast.literal_eval(text), dtype=float)
    except (ValueError, SyntaxError):
        text = text.strip("[]").replace("\n", " ").replace(",", " ")
        return np.fromstring(text, sep=" ")


def load_pair_arrays(n_obs: int) -> tuple[np.ndarray, np.ndarray]:
    usecols = ["ligand", "receptor", "local_moran", "local_geary"]
    for chunk in pd.read_csv(LR_AUTO_CSV, usecols=usecols, chunksize=500):
        direct = chunk[(chunk["ligand"] == LIGAND) & (chunk["receptor"] == RECEPTOR)]
        reverse = chunk[(chunk["ligand"] == RECEPTOR) & (chunk["receptor"] == LIGAND)]
        if not direct.empty:
            row = direct.iloc[0]
        elif not reverse.empty:
            row = reverse.iloc[0]
        else:
            continue

        moran = parse_numeric_list(row["local_moran"])
        geary = parse_numeric_list(row["local_geary"])
        if len(moran) != n_obs or len(geary) != n_obs:
            raise ValueError(
                f"{row['ligand']}_{row['receptor']} local arrays are not aligned: "
                f"moran={len(moran)}, geary={len(geary)}, adata={n_obs}"
            )
        return moran, geary

    raise ValueError(f"Cannot find {SOURCE_PAIR} or {RECEPTOR}_{LIGAND} in {LR_AUTO_CSV}")


def save_spatial(adata, color: str, stem: str, title: str) -> None:
    sc.pl.spatial(
        adata,
        color=color,
        title=title,
        size=1.5,
        cmap="RdYlGn_r",
        alpha_img=0.5,
        show=False,
    )
    plt.savefig(SVG_DIR / f"{stem}.svg", dpi=300, bbox_inches="tight")
    plt.close("all")


def add_state_columns(
    adata,
    prefix: str,
    moran_col: str,
    geary_col: str,
    label: str,
) -> list[dict[str, object]]:
    moran = pd.to_numeric(adata.obs[moran_col], errors="coerce").fillna(0).to_numpy()
    geary = pd.to_numeric(adata.obs[geary_col], errors="coerce").fillna(0).to_numpy()

    rows = []
    for state_key, state_label, mask_fn in STATE_SPECS:
        col = f"{prefix}_{state_key}"
        mask = mask_fn(moran, geary)
        adata.obs[col] = mask.astype(int)
        rows.append(
            {
                "target": label,
                "state": state_label,
                "column": col,
                "spot_count": int(mask.sum()),
                "spot_fraction": float(mask.mean()),
            }
        )
    return rows


def add_log_geary(adata, source_col: str, output_col: str) -> None:
    values = pd.to_numeric(adata.obs[source_col], errors="coerce").fillna(0).to_numpy()
    adata.obs[output_col] = np.log10(values + 1)


def main() -> None:
    ensure_dirs()
    plt.rcParams.update(RC_PARAMS)

    adata = sc.read_h5ad(ADATA_PATH)
    required = [
        f"{RECEPTOR}_moran",
        f"{RECEPTOR}_gearys",
        f"{LIGAND}_moran",
        f"{LIGAND}_gearys",
    ]
    missing = [col for col in required if col not in adata.obs.columns]
    if missing:
        raise ValueError(f"{ADATA_PATH} is missing required obs columns: {missing}")

    pair_moran, pair_geary = load_pair_arrays(adata.n_obs)
    pair_moran_col = f"{DISPLAY_PAIR}_pair_moran"
    pair_geary_col = f"{DISPLAY_PAIR}_pair_gearys"
    pair_geary_log_col = f"{DISPLAY_PAIR}_pair_gearys_log"
    adata.obs[pair_moran_col] = np.nan_to_num(pair_moran, nan=0.0)
    adata.obs[pair_geary_col] = np.nan_to_num(pair_geary, nan=0.0)
    add_log_geary(adata, pair_geary_col, pair_geary_log_col)

    add_log_geary(adata, f"{RECEPTOR}_gearys", f"{RECEPTOR}_gearys_log_redraw")
    add_log_geary(adata, f"{LIGAND}_gearys", f"{LIGAND}_gearys_log_redraw")

    summary_rows = []
    summary_rows.extend(
        add_state_columns(
            adata,
            RECEPTOR,
            f"{RECEPTOR}_moran",
            f"{RECEPTOR}_gearys",
            RECEPTOR,
        )
    )
    summary_rows.extend(
        add_state_columns(
            adata,
            LIGAND,
            f"{LIGAND}_moran",
            f"{LIGAND}_gearys",
            LIGAND,
        )
    )
    summary_rows.extend(
        add_state_columns(
            adata,
            f"{DISPLAY_PAIR}_pair",
            pair_moran_col,
            pair_geary_col,
            f"{DISPLAY_PAIR} pair",
        )
    )
    pd.DataFrame(summary_rows).to_csv(TABLE_DIR / "state_spot_counts.csv", index=False)

    plot_specs = []
    for gene in [RECEPTOR, LIGAND]:
        plot_specs.extend(
            [
                (f"{gene}_moran", f"{gene}_local_moran", f"{gene} local Moran's I"),
                (
                    f"{gene}_gearys_log_redraw",
                    f"{gene}_local_gearys_log",
                    f"{gene} local Geary's C (log10 + 1)",
                ),
            ]
        )
        for state_key, state_label, _ in STATE_SPECS:
            plot_specs.append(
                (
                    f"{gene}_{state_key}",
                    f"{gene}_{state_key}",
                    f"{gene} {state_label}",
                )
            )

    plot_specs.extend(
        [
            (
                pair_moran_col,
                f"{DISPLAY_PAIR}_pair_local_moran",
                f"{DISPLAY_PAIR} pair local Moran's I",
            ),
            (
                pair_geary_log_col,
                f"{DISPLAY_PAIR}_pair_local_gearys_log",
                f"{DISPLAY_PAIR} pair local Geary's C (log10 + 1)",
            ),
        ]
    )
    for state_key, state_label, _ in STATE_SPECS:
        plot_specs.append(
            (
                f"{DISPLAY_PAIR}_pair_{state_key}",
                f"{DISPLAY_PAIR}_pair_{state_key}",
                f"{DISPLAY_PAIR} pair {state_label}",
            )
        )

    for color, stem, title in plot_specs:
        save_spatial(adata, color, stem, title)

    (OUT_DIR / "run_config.json").write_text(
        json.dumps(
            {
                "adata": str(ADATA_PATH),
                "lr_autocorr": str(LR_AUTO_CSV),
                "source_pair": SOURCE_PAIR,
                "display_pair": DISPLAY_PAIR,
                "state_thresholds": {
                    "homogeneous_gradient": "local_moran > 0 and 0 < local_geary <= 0.5",
                    "random_like_normal": "local_moran > 0 and 0.5 < local_geary < 1.5",
                    "heterogeneous_edge": "local_moran > 0 and local_geary > 1.5",
                },
                "geary_display": "log10(local_geary + 1), matching the original notebook-style geary plots",
                "svg_count": len(plot_specs),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print(f"Saved {len(plot_specs)} SVG files to {SVG_DIR}")


if __name__ == "__main__":
    main()
