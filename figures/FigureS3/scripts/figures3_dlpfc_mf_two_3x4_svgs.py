# Submission-aligned copy generated 2026-08-18
# Submission figure: FigureS3
# Role: Figure S3 two separate 3x4 SVG layouts

from __future__ import annotations

import html
import re
import os
from pathlib import Path


SOURCE_DIR = Path(os.environ.get("SPATIALCCCBENCH_SOURCE_PANEL_DIR", Path(__file__).resolve().parents[1] / "results" / "generated" / "individual_panels"))
OUT_DIR = Path(os.environ.get("SPATIALCCCBENCH_OUTPUT_DIR", Path(__file__).resolve().parents[1] / "results" / "generated")) / "composites"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FONT = "Arial"
PAGE_W = 590.0
PAGE_H = 1012.0
MARGIN_X = 18.0
MARGIN_TOP = 18.0
MARGIN_BOTTOM = 18.0
COL_GAP = 8.0
ROW_GAP = 8.0
HEADER_H = 26.0
LR_HEADER_H = 16.0
BOX_PAD_X = 3.0
CELL_W = 174.0
CELL_H = 222.0
BOX_W = CELL_W + 2 * BOX_PAD_X


SVG_ROOT_RE = re.compile(r"<svg\b([^>]*)>", re.IGNORECASE | re.DOTALL)
VIEWBOX_RE = re.compile(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"')
ID_RE = re.compile(r'\bid="([^"]+)"')


DATASETS = [
    {
        "panel": "a",
        "name": "DLPFC",
        "folder": SOURCE_DIR / "DLPFC",
        "lr_label": "SPP1-ITGB1",
        "out_name": "DLPFC_SPP1_ITGB1_12_panels_3x4.svg",
        "headers": [
            "Receptor: ITGB1\nexpression / Geary / Moran",
            "Ligand: SPP1\nexpression / Geary / Moran",
            "SPP1-ITGB1\nper-gene overlap",
        ],
        "rows": [
            [
                "01_DLPFC_SPP1_ITGB1_receptor_raw_expression.svg",
                "02_DLPFC_SPP1_ITGB1_ligand_raw_expression.svg",
                "03_DLPFC_SPP1_ITGB1_per_gene_random_like_normal.svg",
            ],
            [
                "06_DLPFC_SPP1_ITGB1_receptor_local_geary_log.svg",
                "05_DLPFC_SPP1_ITGB1_ligand_local_geary_log.svg",
                "07_DLPFC_SPP1_ITGB1_per_gene_homogeneous_gradient.svg",
            ],
            [
                "10_DLPFC_SPP1_ITGB1_receptor_local_moran.svg",
                "09_DLPFC_SPP1_ITGB1_ligand_local_moran.svg",
                "11_DLPFC_SPP1_ITGB1_per_gene_heterogeneous_edge.svg",
            ],
            [
                "04_DLPFC_SPP1_ITGB1_lr_integrated_lr_geary_0p5_1p5.svg",
                "08_DLPFC_SPP1_ITGB1_lr_integrated_lr_geary_0_0p5.svg",
                "12_DLPFC_SPP1_ITGB1_lr_integrated_lr_geary_1p5_inf.svg",
            ],
        ],
    },
    {
        "panel": "b",
        "name": "MF",
        "folder": SOURCE_DIR / "MF",
        "lr_label": "Pnoc-Oprl1",
        "out_name": "MF_Pnoc_Oprl1_12_panels_3x4.svg",
        "headers": [
            "Receptor: Oprl1\nexpression / Geary / Moran",
            "Ligand: Pnoc\nexpression / Geary / Moran",
            "Pnoc-Oprl1\nper-gene overlap",
        ],
        "rows": [
            [
                "13_MF_Pnoc_Oprl1_receptor_raw_expression.svg",
                "14_MF_Pnoc_Oprl1_ligand_raw_expression.svg",
                "15_MF_Pnoc_Oprl1_per_gene_random_like_normal.svg",
            ],
            [
                "18_MF_Pnoc_Oprl1_receptor_local_geary_log.svg",
                "17_MF_Pnoc_Oprl1_ligand_local_geary_log.svg",
                "19_MF_Pnoc_Oprl1_per_gene_homogeneous_gradient.svg",
            ],
            [
                "22_MF_Pnoc_Oprl1_receptor_local_moran.svg",
                "21_MF_Pnoc_Oprl1_ligand_local_moran.svg",
                "23_MF_Pnoc_Oprl1_per_gene_heterogeneous_edge.svg",
            ],
            [
                "16_MF_Pnoc_Oprl1_lr_integrated_lr_geary_0p5_1p5.svg",
                "20_MF_Pnoc_Oprl1_lr_integrated_lr_geary_0_0p5.svg",
                "24_MF_Pnoc_Oprl1_lr_integrated_lr_geary_1p5_inf.svg",
            ],
        ],
    },
]


def read_svg(path: Path, prefix: str) -> tuple[str, float, float]:
    text = path.read_text(encoding="utf-8")
    match = SVG_ROOT_RE.search(text)
    if match is None:
        raise ValueError(f"Missing svg root: {path}")
    viewbox = VIEWBOX_RE.search(match.group(0))
    if viewbox is None:
        raise ValueError(f"Missing viewBox: {path}")
    width = float(viewbox.group(1))
    height = float(viewbox.group(2))
    inner = text[match.end() : text.rfind("</svg>")]
    inner = re.sub(r"<metadata>.*?</metadata>", "", inner, flags=re.DOTALL)
    ids = set(ID_RE.findall(inner))
    if ids:
        id_map = {old_id: f"{prefix}_{old_id}" for old_id in ids}
        inner = re.sub(r'\bid="([^"]+)"', lambda m: f'id="{id_map.get(m.group(1), m.group(1))}"', inner)
        inner = re.sub(r"url\(#([^)]+)\)", lambda m: f"url(#{id_map.get(m.group(1), m.group(1))})", inner)
        inner = re.sub(
            r'(xlink:href|href)="#([^"]+)"',
            lambda m: f'{m.group(1)}="#{id_map.get(m.group(2), m.group(2))}"',
            inner,
        )
    return inner, width, height


def multiline_text(label: str, x: float, y: float, size: float, line_gap: float = 8.0) -> str:
    lines = label.split("\n")
    start_y = y - (len(lines) - 1) * line_gap / 2
    return "\n".join(
        [
            f'<text x="{x:.3f}" y="{start_y + idx * line_gap:.3f}" '
            f'font-family="{FONT}" font-size="{size}" text-anchor="middle" '
            f'dominant-baseline="middle">{html.escape(line)}</text>'
            for idx, line in enumerate(lines)
        ]
    )


def panel_position(row: int, col: int) -> tuple[float, float]:
    x = MARGIN_X + col * (BOX_W + COL_GAP) + BOX_PAD_X
    y = MARGIN_TOP + HEADER_H + row * (CELL_H + ROW_GAP)
    if row == 3:
        y += LR_HEADER_H
    return x, y


def build_dataset_svg(config: dict[str, object]) -> Path:
    parts: list[str] = [f'<rect x="0" y="0" width="{PAGE_W:.3f}" height="{PAGE_H:.3f}" fill="#FFFFFF"/>']

    top_y = MARGIN_TOP
    box_h = HEADER_H + 3 * CELL_H + 2 * ROW_GAP
    bottom_y = top_y + box_h + ROW_GAP
    bottom_h = LR_HEADER_H + CELL_H

    for col, label in enumerate(config["headers"]):
        bx = MARGIN_X + col * (BOX_W + COL_GAP)
        parts.append(
            f'<rect x="{bx:.3f}" y="{top_y:.3f}" width="{BOX_W:.3f}" height="{box_h:.3f}" '
            f'fill="none" stroke="#4D4D4D" stroke-width="0.55"/>'
        )
        parts.append(
            f'<rect x="{bx:.3f}" y="{top_y:.3f}" width="{BOX_W:.3f}" height="{HEADER_H:.3f}" '
            f'fill="#F2F2F2" stroke="#4D4D4D" stroke-width="0.55"/>'
        )
        parts.append(multiline_text(str(label), bx + BOX_W / 2, top_y + HEADER_H / 2, 6.4))

    bottom_x = MARGIN_X
    bottom_w = 3 * BOX_W + 2 * COL_GAP
    parts.append(
        f'<rect x="{bottom_x:.3f}" y="{bottom_y:.3f}" width="{bottom_w:.3f}" height="{bottom_h:.3f}" '
        f'fill="none" stroke="#4D4D4D" stroke-width="0.55"/>'
    )
    parts.append(
        f'<rect x="{bottom_x:.3f}" y="{bottom_y:.3f}" width="{bottom_w:.3f}" height="{LR_HEADER_H:.3f}" '
        f'fill="#F2F2F2" stroke="#4D4D4D" stroke-width="0.55"/>'
    )
    parts.append(
        f'<text x="{bottom_x + bottom_w / 2:.3f}" y="{bottom_y + LR_HEADER_H / 2:.3f}" '
        f'font-family="{FONT}" font-size="6.8" text-anchor="middle" dominant-baseline="middle">'
        f'{html.escape(str(config["lr_label"]))} LR-integrated Geary</text>'
    )

    rows = config["rows"]
    folder = Path(config["folder"])
    for row_idx, row_files in enumerate(rows):
        for col_idx, filename in enumerate(row_files):
            svg, width, height = read_svg(folder / str(filename), f"{config['name']}_{row_idx + 1}_{col_idx + 1}")
            x, y = panel_position(row_idx, col_idx)
            scale = min(CELL_W / width, CELL_H / height)
            tx = x + (CELL_W - width * scale) / 2
            ty = y + (CELL_H - height * scale) / 2
            parts.append(f'<g transform="translate({tx:.3f},{ty:.3f}) scale({scale:.6f})">\n{svg}\n</g>')

    svg_text = (
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{PAGE_W:.3f}pt" height="{PAGE_H:.3f}pt" viewBox="0 0 {PAGE_W:.3f} {PAGE_H:.3f}">\n'
        f'{"".join(parts)}\n</svg>\n'
    )
    out_path = OUT_DIR / str(config["out_name"])
    out_path.write_text(svg_text, encoding="utf-8")
    return out_path


def main() -> None:
    for config in DATASETS:
        out_path = build_dataset_svg(config)
        print(out_path)


if __name__ == "__main__":
    main()
