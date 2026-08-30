from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


EXPECTED = ["Figure1", "Figure2", "Figure3", "Figure4", "Figure5"] + [f"FigureS{i}" for i in range(1, 10)]
COMPUTATIONAL = {"Figure2", "Figure3", "Figure4", "Figure5", "FigureS3"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def indexed_csv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}
        index_column = reader.fieldnames[0]
        return {row[index_column]: row for row in reader}


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    errors = []
    for name in EXPECTED:
        folder = root / "figures" / name
        if not folder.is_dir():
            errors.append(f"missing directory: {folder}")
            continue
        if len(list(folder.glob("run.py"))) != 1:
            errors.append(f"missing run.py: {name}")
        notebooks = list(folder.glob("*_demo.ipynb"))
        if len(notebooks) != 1:
            errors.append(f"missing demo notebook: {name}")
        else:
            try:
                notebook = json.loads(notebooks[0].read_text(encoding="utf-8"))
                code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
                if any(cell.get("outputs") or cell.get("execution_count") is not None for cell in code_cells):
                    errors.append(f"notebook contains stale execution output: {name}")
                source = "".join("".join(cell.get("source", [])) for cell in code_cells)
                if name in COMPUTATIONAL and "--reproduce" not in source:
                    errors.append(f"notebook does not call the reproducible shortcut: {name}")
            except Exception as exc:
                errors.append(f"invalid demo notebook {notebooks[0]}: {exc}")
        result = folder / "results" / "published" / f"{name}.svg"
        if not result.exists() or result.stat().st_size == 0:
            errors.append(f"missing published result: {result}")
        config = folder / "figure_config.json"
        try:
            json.loads(config.read_text(encoding="utf-8"))
        except Exception as exc:
            errors.append(f"invalid config {config}: {exc}")
        if name in COMPUTATIONAL:
            generated = folder / "results" / "generated"
            if not any(generated.rglob("*.svg")):
                errors.append(f"missing script-generated SVG results: {name}")
            manifest_path = folder / "results" / "result_manifest.json"
            try:
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                if not manifest.get("redraw_requested") or not manifest.get("generated_files"):
                    errors.append(f"result manifest does not describe a completed redraw: {name}")
            except Exception as exc:
                errors.append(f"invalid result manifest {manifest_path}: {exc}")

    counts_path = root / "figures" / "Figure3" / "inputs" / "boundary_ligand_counts.csv"
    counts_metadata_path = counts_path.with_suffix(".json")
    try:
        counts_metadata = json.loads(counts_metadata_path.read_text(encoding="utf-8"))
        if sha256(counts_path) != counts_metadata["output_sha256"]:
            errors.append("Figure3 compact-count hash does not match its metadata")
    except Exception as exc:
        errors.append(f"invalid Figure3 compact-count provenance: {exc}")

    figure4_summary = root / "figures" / "Figure4" / "inputs" / "result_summary.csv"
    figure5_summary = root / "figures" / "Figure5" / "inputs" / "result_summary.csv"
    if figure4_summary.read_bytes() != figure5_summary.read_bytes():
        errors.append("Figure4 and Figure5 result_summary.csv files are not synchronized")

    try:
        summary = indexed_csv(figure5_summary)
        recomputed = indexed_csv(
            root / "figures" / "Figure5" / "inputs" / "boundary_metrics_recomputed.csv"
        )
        boundary_columns = ("dropdown_edge", "equal_edge", "gradient_edge")
        for tool, values in recomputed.items():
            for column in boundary_columns:
                if not math.isclose(
                    float(values[column]),
                    float(summary[tool][column]),
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ):
                    errors.append(f"Figure5 boundary mismatch: {tool}/{column}")
    except Exception as exc:
        errors.append(f"invalid Figure5 boundary audit: {exc}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"Validated {len(EXPECTED)} Figure directories.")


if __name__ == "__main__":
    main()
