from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract CCC results and calculate precision/recall/F1.")
    parser.add_argument("--dataset", default="st_lymphnode")
    parser.add_argument("--spot-info", type=Path, required=True)
    parser.add_argument("--result-root", type=Path, default=Path("result"))
    parser.add_argument("--tools", nargs="+", default=["CellphoneDB", "Baseline_1", "Baseline_2"])
    parser.add_argument("--ground-truth-level", type=int, default=3)
    parser.add_argument("--output", type=Path, default=Path("tutorial_outputs/accuracy_metrics.csv"))
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    os.chdir(repo_root)
    from toolkit.Cal_precision import (
        get_LR_from_per_cellpair_precision,
        obtain_LR_cell_list,
        obtain_all_result,
    )
    from toolkit.preprocess import extract_result

    if args.result_root.name != "result":
        raise ValueError("--result-root must be the repository result directory because toolkit paths are relative.")
    extracted = extract_result(args.tools, args.dataset, spot_info=str(args.spot_info))
    lr_cell_list = obtain_LR_cell_list(extracted)
    combined = obtain_all_result(lr_cell_list)
    metrics = get_LR_from_per_cellpair_precision(combined, ground_truth_level=args.ground_truth_level)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.output, index=False)
    print(f"Wrote {len(metrics)} metric rows to {args.output}")


if __name__ == "__main__":
    main()
