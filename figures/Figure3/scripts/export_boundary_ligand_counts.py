from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import pandas as pd

from figure3_abc_boundary_gradient_equalization_signal_loss import LEVEL_LIST, TOOL_LIST
from figure3_preprocess_helper import extract_result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ligand_counts(frame: pd.DataFrame) -> pd.Series:
    pairs = frame["LR_pairs"].astype(str).str.replace("HLA_", "HLA-", regex=False)
    return pairs.str.split("_", n=1, expand=True)[0].value_counts()


def main() -> None:
    parser = argparse.ArgumentParser(description="Export compact Figure3 ligand counts.")
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data_root = args.data_root.resolve()
    previous = Path.cwd()
    rows = []
    try:
        os.chdir(data_root)
        for level in LEVEL_LIST:
            print(f"Extracting {level}...", flush=True)
            extracted = extract_result(
                TOOL_LIST,
                level,
                spot_info="./dataset/DLPFC_cell_info.csv",
            )
            for tool, frame in extracted.items():
                for ligand, count in ligand_counts(frame).items():
                    rows.append(
                        {
                            "scenario": level,
                            "tool": tool,
                            "ligand": ligand,
                            "count": int(count),
                        }
                    )
    finally:
        os.chdir(previous)

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows).sort_values(["scenario", "tool", "ligand"])
    result.to_csv(output, index=False)
    metadata = {
        "generated_from": "external-local-data",
        "rows": len(result),
        "scenarios": LEVEL_LIST,
        "tools": TOOL_LIST,
        "output_sha256": sha256(output),
    }
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Saved {len(result)} rows to {output}")


if __name__ == "__main__":
    main()
