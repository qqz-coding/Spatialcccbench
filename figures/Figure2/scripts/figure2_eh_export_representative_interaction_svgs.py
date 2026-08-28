# Submission-aligned copy generated 2026-08-18
# Submission figure: Figure2
# Role: Figure 2E-H export individual interaction SVGs

from __future__ import annotations

import ast
import json
import shutil
import warnings
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


warnings.filterwarnings("ignore")

REPO_ROOT = Path(os.environ.get("SPATIALCCCBENCH_REPO_ROOT", Path(__file__).resolve().parents[3]))
PROJECT_DIR = Path(os.environ.get("SPATIALCCCBENCH_DATA_ROOT", REPO_ROOT))
OUT_DIR = Path(os.environ.get("SPATIALCCCBENCH_OUTPUT_DIR", Path(__file__).resolve().parents[1] / "results" / "generated"))
SVG_DIR = OUT_DIR / "svg"
TABLE_DIR = OUT_DIR / "tables"

ADATA_PATH = PROJECT_DIR / "dataset" / "st_lymphnode.h5ad"
LIGAND_RESULT = PROJECT_DIR / "result_l_1_full.csv"
RECEPTOR_RESULT = PROJECT_DIR / "result_r_1_full.csv"

EXAMPLES = [
    {
        "label": "local_diffusion",
        "category": "heterogeneous_edge",
        "ligand": "IGF1",
        "receptor": "ITGB4",
        "mask_name": "upheaval",
        "interaction_prefix": "upheaval",
        "filename": "01_local_diffusion_IGF1_ITGB4_interaction.svg",
    },
    {
        "label": "local_colocalization",
        "category": "homogeneous_gradient",
        "ligand": "CXCL10",
        "receptor": "CXCR3",
        "mask_name": "colocalization",
        "interaction_prefix": "colocalization",
        "filename": "02_local_colocalization_CXCL10_CXCR3_interaction.svg",
    },
    {
        "label": "global_colocalization",
        "category": "random_like_normal",
        "ligand": "PLAU",
        "receptor": "PLAUR",
        "mask_name": "colocalization",
        "interaction_prefix": "colocalization",
        "filename": "03_global_colocalization_PLAU_PLAUR_interaction.svg",
    },
    {
        "label": "weak_edge_no_interaction",
        "category": "edge_no_interaction",
        "ligand": "SLC1A6",
        "receptor": "GRIA1",
        "mask_name": "upheaval",
        "interaction_prefix": "upheaval",
        "filename": "04_weak_edge_no_interaction_SLC1A6_GRIA1_interaction.svg",
    },
]


def reset_output_dirs() -> None:
    for path in [SVG_DIR, TABLE_DIR]:
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True, exist_ok=True)


def clean_symbol(value: object) -> str:
    return str(value).strip().strip("'\"")


def parse_array(value: object) -> np.ndarray:
    return np.asarray(ast.literal_eval(str(value)), dtype=float)


def load_gene_results(path: Path, gene_col: str, genes: set[str]) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=[gene_col, "local_morani", "local_gearysC"])
    df[gene_col] = df[gene_col].map(clean_symbol)
    df = df[df[gene_col].isin(genes)].copy()
    missing = sorted(genes - set(df[gene_col]))
    if missing:
        raise ValueError(f"Missing {gene_col} rows in {path}: {missing}")
    df["local_morani"] = df["local_morani"].map(parse_array)
    df["local_gearysC"] = df["local_gearysC"].map(parse_array)
    return df.set_index(gene_col)


def mask_for_category(local_moran: np.ndarray, local_geary: np.ndarray, category: str) -> np.ndarray:
    if category == "heterogeneous_edge":
        return (local_moran >= 0) & (local_geary >= 1.5)
    if category == "homogeneous_gradient":
        return (local_moran >= 0) & (local_geary > 0) & (local_geary < 0.5)
    if category == "random_like_normal":
        return (local_moran >= 0) & (local_geary > 0.5) & (local_geary < 1.5)
    if category == "edge_no_interaction":
        return (local_moran > 0) & (local_geary >= 1.5)
    raise ValueError(category)


def add_mask_columns(
    adata: sc.AnnData,
    gene: str,
    result_df: pd.DataFrame,
    category: str,
    mask_name: str,
) -> np.ndarray:
    local_moran = np.asarray(result_df.loc[gene, "local_morani"], dtype=float)
    local_geary = np.asarray(result_df.loc[gene, "local_gearysC"], dtype=float)
    if local_moran.shape[0] != adata.n_obs or local_geary.shape[0] != adata.n_obs:
        raise ValueError(f"{gene} local array length does not match adata.n_obs")
    adata.obs[f"{gene}_moran"] = local_moran
    adata.obs[f"{gene}_gearys"] = local_geary
    adata.obs[f"{gene}_gearys(log)"] = np.log10(np.clip(local_geary, a_min=0, a_max=None) + 1)
    mask = mask_for_category(local_moran, local_geary, category)
    adata.obs[f"{gene}_{mask_name}"] = 0
    adata.obs.loc[mask, f"{gene}_{mask_name}"] = 1
    return mask


def add_interaction_column(
    adata: sc.AnnData,
    ligand: str,
    receptor: str,
    ligand_mask: np.ndarray,
    receptor_mask: np.ndarray,
    interaction_prefix: str,
) -> str:
    column = f"{ligand}_{receptor}_interaction"
    adata.obs[column] = "other"
    adata.obs.loc[ligand_mask, column] = f"ligand_{interaction_prefix}"
    adata.obs.loc[receptor_mask, column] = f"receptor_{interaction_prefix}"
    adata.obs.loc[ligand_mask & receptor_mask, column] = "interaction"
    return column


def save_notebook_style_svg(adata: sc.AnnData, color: str, out_path: Path) -> None:
    # Keep the plotting call aligned with the saved cells in Figure1-A.ipynb.
    sc.pl.spatial(
        adata,
        color=color,
        size=1.5,
        cmap="RdYlGn_r",
        alpha_img=0.5,
        show=False,
    )
    plt.savefig(out_path, dpi=300, bbox_inches="tight", format="svg")
    plt.close()


def run() -> None:
    reset_output_dirs()

    ligand_genes = {example["ligand"] for example in EXAMPLES}
    receptor_genes = {example["receptor"] for example in EXAMPLES}
    ligand_results = load_gene_results(LIGAND_RESULT, "ligand", ligand_genes)
    receptor_results = load_gene_results(RECEPTOR_RESULT, "receptor", receptor_genes)
    adata = sc.read_h5ad(ADATA_PATH)

    records: list[dict[str, object]] = []
    for example in EXAMPLES:
        ligand = str(example["ligand"])
        receptor = str(example["receptor"])
        ligand_mask = add_mask_columns(
            adata,
            ligand,
            ligand_results,
            str(example["category"]),
            str(example["mask_name"]),
        )
        receptor_mask = add_mask_columns(
            adata,
            receptor,
            receptor_results,
            str(example["category"]),
            str(example["mask_name"]),
        )
        interaction_col = add_interaction_column(
            adata,
            ligand,
            receptor,
            ligand_mask,
            receptor_mask,
            str(example["interaction_prefix"]),
        )

        out_path = SVG_DIR / str(example["filename"])
        save_notebook_style_svg(adata, interaction_col, out_path)

        records.append(
            {
                "label": example["label"],
                "LR": f"{ligand}_{receptor}",
                "category": example["category"],
                "ligand_mask_count": int(ligand_mask.sum()),
                "receptor_mask_count": int(receptor_mask.sum()),
                "interaction_count": int((ligand_mask & receptor_mask).sum()),
                "svg": str(out_path),
            }
        )

    pd.DataFrame(records).to_csv(TABLE_DIR / "four_svg_summary.csv", index=False)
    (OUT_DIR / "run_config.json").write_text(
        json.dumps(
            {
                "adata": str(ADATA_PATH),
                "ligand_result": str(LIGAND_RESULT),
                "receptor_result": str(RECEPTOR_RESULT),
                "array_parser": "ast.literal_eval",
                "plotting_call": "sc.pl.spatial(adata, color=..., size=1.5, cmap='RdYlGn_r', alpha_img=0.5, show=False); plt.savefig(..., format='svg')",
                "examples": EXAMPLES,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    run()
