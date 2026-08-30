from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    rows: list[dict[str, object]] = []
    for figure_dir in sorted((REPO_ROOT / "figures").glob("Figure*")):
        config_path = figure_dir / "figure_config.json"
        if not config_path.exists():
            continue
        config = json.loads(config_path.read_text(encoding="utf-8"))
        inputs_dir = figure_dir / "inputs"
        if inputs_dir.exists():
            for path in sorted(inputs_dir.rglob("*")):
                if not path.is_file() or path.suffix.lower() not in {".csv", ".json", ".h5ad"}:
                    continue
                rows.append(
                    {
                        "figure": config["figure"],
                        "status": config["status"],
                        "result_kind": "input",
                        "path": path.relative_to(REPO_ROOT).as_posix(),
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )
        for result_kind in ("published", "generated"):
            result_dir = figure_dir / "results" / result_kind
            if not result_dir.exists():
                continue
            manifest_excludes = (
                set(config.get("generated_manifest_exclude", []))
                if result_kind == "generated"
                else set()
            )
            for path in sorted(result_dir.rglob("*")):
                if not path.is_file() or path.suffix.lower() not in {".svg", ".csv", ".json"}:
                    continue
                if manifest_excludes.intersection(path.relative_to(result_dir).parts):
                    continue
                rows.append(
                    {
                        "figure": config["figure"],
                        "status": config["status"],
                        "result_kind": result_kind,
                        "path": path.relative_to(REPO_ROOT).as_posix(),
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )

    output = REPO_ROOT / "results_manifest.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["figure", "status", "result_kind", "path", "bytes", "sha256"],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} result records to {output}")


if __name__ == "__main__":
    main()
