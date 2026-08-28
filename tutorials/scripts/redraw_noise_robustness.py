from __future__ import annotations

import subprocess
import sys
from pathlib import Path


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    subprocess.run([sys.executable, str(root / "figures" / "Figure4" / "run.py"), "--reproduce"], cwd=root, check=True)
