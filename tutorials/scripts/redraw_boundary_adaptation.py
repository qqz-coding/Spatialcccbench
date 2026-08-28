from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    subprocess.run([sys.executable, str(root / "figures" / "Figure3" / "run.py"), "--reproduce", "--data-root", str(args.data_root)], cwd=root, check=True)
