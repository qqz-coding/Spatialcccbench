from __future__ import annotations

import argparse
import shutil
from datetime import datetime
import os
from pathlib import Path

import pandas as pd

from plotting_module import NOISE_ORDER, plot_all_noise_resistance


ROOT = Path(os.environ.get("SPATIALCCCBENCH_DATA_ROOT", Path(__file__).resolve().parents[3]))
DEFAULT_NOISE = ROOT / "processed_noise.csv"
DEFAULT_SUMMARY = ROOT / "result_summary.csv"
DEFAULT_OUT = ROOT / "figure4_noise_resistance_processed_20260827"


def read_noise(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    first = df.columns[0]
    if first.startswith("Unnamed") or first == "":
        df = df.rename(columns={first: "noise"})
    required = {"noise", "precision", "recall", "F1", "tool"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    df["noise"] = df["noise"].astype(str)
    df["tool"] = df["tool"].astype(str)
    return df


def update_result_summary(noise_df: pd.DataFrame, summary_path: Path, out_dir: Path) -> Path:
    summary = pd.read_csv(summary_path, index_col=0)
    summary.index = summary.index.astype(str)
    f1 = noise_df.pivot(index="tool", columns="noise", values="F1").reindex(columns=NOISE_ORDER)
    f1.index = f1.index.str.replace("stLearn_without_spotmixture", "stLearn*", regex=False)
    missing_tools = sorted(set(summary.index) - set(f1.index))
    if missing_tools:
        raise ValueError(f"Noise table is missing summary tools: {missing_tools}")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = out_dir / "backups"
    backup_dir.mkdir(parents=True, exist_ok=True)
    backup_path = backup_dir / f"result_summary.before_processed_noise_{timestamp}.csv"
    shutil.copy2(summary_path, backup_path)

    changes = []
    for noise in NOISE_ORDER:
        column = f"{noise}_F1"
        old = summary[column].astype(float).copy() if column in summary else pd.Series(index=summary.index, dtype=float)
        new = f1[noise].reindex(summary.index).astype(float)
        summary[column] = new
        changes.append(
            pd.DataFrame(
                {
                    "tool": summary.index,
                    "noise": noise,
                    "old_F1": old.reindex(summary.index).to_numpy(),
                    "new_F1": new.to_numpy(),
                    "delta": (new - old.reindex(summary.index)).to_numpy(),
                }
            )
        )

    summary.to_csv(summary_path)
    summary.to_csv(out_dir / "result_summary.updated.csv")
    pd.concat(changes, ignore_index=True).to_csv(out_dir / "result_summary_noise_update.csv", index=False)
    return backup_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Render processed noise metrics as individual radar SVGs.")
    parser.add_argument("--noise-csv", type=Path, default=DEFAULT_NOISE)
    parser.add_argument("--summary-csv", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--skip-summary-update", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    noise_df = read_noise(args.noise_csv)
    noise_df.to_csv(args.out_dir / "processed_noise.used.csv", index=False)
    paths, legend = plot_all_noise_resistance(noise_df, args.out_dir / "svg")
    print(f"Wrote {len(paths)} tool SVGs and legend: {legend}")

    if not args.skip_summary_update:
        backup = update_result_summary(noise_df, args.summary_csv, args.out_dir)
        print(f"Updated {args.summary_csv}; backup: {backup}")


if __name__ == "__main__":
    main()
