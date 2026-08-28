from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Redraw spatial figures from precomputed h5ad/CSV results.")
    parser.add_argument("figure", choices=["Figure2", "FigureS3"])
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--step", action="append", type=int)
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    command = [sys.executable, str(repo_root / "figures" / args.figure / "run.py"), "--reproduce", "--data-root", str(args.data_root)]
    for step in args.step or []:
        command.extend(["--step", str(step)])
    subprocess.run(command, cwd=repo_root, check=True)


if __name__ == "__main__":
    main()
