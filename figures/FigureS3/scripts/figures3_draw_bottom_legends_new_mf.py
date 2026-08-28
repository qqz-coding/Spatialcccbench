# Submission-aligned copy generated 2026-08-18
# Submission figure: FigureS3
# Role: Intermediate Figure S3 bottom legends layout with updated MF data

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import scanpy as sc


BASE_SCRIPT = Path(__file__).with_name("figures3_draw_compact_layout.py")
spec = importlib.util.spec_from_file_location("compact_layout", BASE_SCRIPT)
compact = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(compact)


OUT_DIR = Path(os.environ.get("SPATIALCCCBENCH_OUTPUT_DIR", Path(__file__).resolve().parents[1] / "results" / "generated"))
REPO_ROOT = Path(os.environ.get("SPATIALCCCBENCH_REPO_ROOT", Path(__file__).resolve().parents[3]))
PROJECT_DIR = Path(os.environ.get("SPATIALCCCBENCH_DATA_ROOT", REPO_ROOT))
compact.OUT_DIR = OUT_DIR
compact.CACHE_DIR = OUT_DIR / "h5ad_readable_cache"
compact.PANEL_DIR = OUT_DIR / "individual_panels"
compact.COMPOSITE_DIR = OUT_DIR / "composites"
compact.TABLE_DIR = OUT_DIR / "tables"

compact.DATASETS = [
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


def add_bottom_legend(ax: plt.Axes, categories: list[str], palette: dict[str, str]) -> None:
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
        loc="upper center",
        bbox_to_anchor=(0.5, -0.155),
        ncol=2,
        frameon=False,
        handletextpad=0.25,
        columnspacing=0.70,
        borderaxespad=0,
        labelspacing=0.18,
        fontsize=5.6,
    )


def plot_continuous_bottom(
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
        fig, ax = plt.subplots(figsize=(2.65, 2.65))
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
        cmap=compact.CONTINUOUS_CMAP,
        colorbar_loc=None,
        use_raw=use_raw,
    )
    mappable = next((coll for coll in reversed(ax.collections) if coll.get_array() is not None), None)
    if mappable is not None:
        cax = ax.inset_axes([0.0, -0.285, 1.0, 0.045], transform=ax.transAxes)
        colorbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        colorbar.ax.tick_params(labelsize=5.8, width=0.4, length=2.0, pad=1.0)
        colorbar.outline.set_linewidth(0.4)
    if created and output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


def plot_categorical_bottom(
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
        fig, ax = plt.subplots(figsize=(2.85, 2.65))
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
    add_bottom_legend(ax, categories, palette)
    if created and output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


def render_spec_bottom(spec: dict[str, object], ax: plt.Axes | None = None, output_path: Path | None = None) -> None:
    if spec["kind"] == "continuous":
        plot_continuous_bottom(
            spec["adata"],
            str(spec["color"]),
            title=str(spec["title"]),
            spot_size=float(spec["spot_size"]),
            use_raw=spec["use_raw"],
            ax=ax,
            output_path=output_path,
        )
    else:
        plot_categorical_bottom(
            spec["adata"],
            str(spec["color"]),
            spec["categories"],
            spec["palette"],
            title=str(spec["title"]),
            spot_size=float(spec["spot_size"]),
            ax=ax,
            output_path=output_path,
        )


def save_dataset_composite_bottom(config: dict[str, object], rows: list[list[dict[str, object]]]) -> None:
    name = str(config["name"])
    lr_pair = str(config["lr_pair"])
    fig = plt.figure(figsize=(12.0, 7.45))
    gs = GridSpec(3, 4, figure=fig, wspace=0.25, hspace=0.54)
    for r in range(3):
        for c in range(4):
            render_spec_bottom(rows[r][c], ax=fig.add_subplot(gs[r, c]))
    fig.subplots_adjust(left=0.045, right=0.975, top=0.965, bottom=0.070)
    fig.savefig(compact.COMPOSITE_DIR / f"{name}_{lr_pair}_12_panels_bottom.svg", dpi=300, bbox_inches="tight")
    fig.savefig(compact.COMPOSITE_DIR / f"{name}_{lr_pair}_12_panels_bottom.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_full_composite_bottom(dataset_rows: list[tuple[dict[str, object], list[list[dict[str, object]]]]]) -> None:
    fig = plt.figure(figsize=(12.2, 14.95))
    gs = GridSpec(
        7,
        4,
        figure=fig,
        height_ratios=[1, 1, 1, 0.10, 1, 1, 1],
        wspace=0.25,
        hspace=0.54,
    )
    row_offsets = [0, 4]
    panel_y = [0.982, 0.492]
    for dataset_idx, (config, rows) in enumerate(dataset_rows):
        offset = row_offsets[dataset_idx]
        fig.text(0.014, panel_y[dataset_idx], str(config["panel"]), fontsize=22, fontfamily="Arial")
        for r in range(3):
            for c in range(4):
                render_spec_bottom(rows[r][c], ax=fig.add_subplot(gs[offset + r, c]))
    fig.subplots_adjust(left=0.045, right=0.975, top=0.982, bottom=0.035)
    fig.savefig(compact.COMPOSITE_DIR / "FigureS3_like_24_bottom_legends_new_mf.svg", dpi=300, bbox_inches="tight")
    fig.savefig(compact.COMPOSITE_DIR / "FigureS3_like_24_bottom_legends_new_mf.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    compact.ensure_dirs()
    plt.rcParams.update(compact.RC_PARAMS)
    compact.render_spec = render_spec_bottom
    compact.save_dataset_composite = save_dataset_composite_bottom
    compact.save_full_composite = save_full_composite_bottom
    compact.main()


if __name__ == "__main__":
    main()
