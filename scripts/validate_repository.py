from __future__ import annotations

import json
from pathlib import Path


EXPECTED = ["Figure1", "Figure2", "Figure3", "Figure4", "Figure5"] + [f"FigureS{i}" for i in range(1, 10)]


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
        if len(list(folder.glob("*_demo.ipynb"))) != 1:
            errors.append(f"missing demo notebook: {name}")
        result = folder / "results" / "published" / f"{name}.svg"
        if not result.exists() or result.stat().st_size == 0:
            errors.append(f"missing published result: {result}")
        config = folder / "figure_config.json"
        try:
            json.loads(config.read_text(encoding="utf-8"))
        except Exception as exc:
            errors.append(f"invalid config {config}: {exc}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"Validated {len(EXPECTED)} Figure directories.")


if __name__ == "__main__":
    main()
