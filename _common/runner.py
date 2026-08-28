from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_figure(figure_dir: Path) -> None:
    config = json.loads((figure_dir / "figure_config.json").read_text(encoding="utf-8"))
    parser = argparse.ArgumentParser(description=config["caption"])
    parser.add_argument("--reproduce", action="store_true", help="Run available redraw steps.")
    parser.add_argument("--step", action="append", type=int, help="Run selected 1-based step(s) only.")
    parser.add_argument("--data-root", type=Path, help="Directory containing external h5ad/result inputs.")
    parser.add_argument("--python", default=sys.executable, help="Python interpreter for redraw scripts.")
    args = parser.parse_args()

    repo_root = figure_dir.parents[1]
    generated_dir = figure_dir / "results" / "generated"
    generated_dir.mkdir(parents=True, exist_ok=True)
    data_root = args.data_root or Path(os.environ.get("SPATIALCCCBENCH_DATA_ROOT", repo_root / "data"))
    env = dict(os.environ)
    env.update(
        {
            "MPLBACKEND": "Agg",
            "PYTHONUNBUFFERED": "1",
            "SPATIALCCCBENCH_REPO_ROOT": str(repo_root),
            "SPATIALCCCBENCH_DATA_ROOT": str(data_root),
            "SPATIALCCCBENCH_OUTPUT_DIR": str(generated_dir),
        }
    )
    summary = figure_dir / "inputs" / "result_summary.csv"
    if summary.exists():
        env["SPATIALCCCBENCH_SUMMARY_CSV"] = str(summary)

    commands = config.get("steps", [])
    if args.reproduce:
        selected = commands
        if args.step:
            invalid = [n for n in args.step if n < 1 or n > len(commands)]
            if invalid:
                raise ValueError(f"Invalid --step values {invalid}; available range is 1-{len(commands)}")
            selected = [commands[n - 1] for n in args.step]
        for step in selected:
            step_args = [str(part).replace("{output}", str(generated_dir)) for part in step.get("args", [])]
            command = [args.python, str(figure_dir / step["script"]), *step_args]
            print("Running:", subprocess.list2cmdline(command), flush=True)
            subprocess.run(command, cwd=data_root if step.get("cwd") == "data-root" else repo_root, env=env, check=True)
    else:
        print(f"{config['figure']}: reference/result files are already prepared; use --reproduce to redraw.")

    published = figure_dir / config["published_result"]
    manifest = {
        "figure": config["figure"],
        "status": config["status"],
        "caption": config["caption"],
        "updated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "published_result": str(published),
        "published_result_sha256": sha256(published) if published.exists() else None,
        "data_root": str(data_root),
        "redraw_requested": args.reproduce,
        "steps": [step["name"] for step in commands],
        "notes": config["notes"],
    }
    (figure_dir / "results" / "result_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"Published result: {published}")
