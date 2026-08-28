from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


NOISE_COLUMNS = [
    "non_specific_F1",
    "up_F1",
    "offset_F1",
    "overlap_F1",
    "dropout_F1",
    "lack_F1",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Update only the six noise columns in Figure5 result_summary.csv."
    )
    parser.add_argument("--base-summary", type=Path, required=True)
    parser.add_argument("--noise-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base = pd.read_csv(args.base_summary).set_index("tool")
    noise = pd.read_csv(args.noise_summary).set_index("tool")
    missing = [column for column in NOISE_COLUMNS if column not in noise.columns]
    if missing:
        raise ValueError(f"Noise summary is missing columns: {missing}")
    if set(base.index) != set(noise.index):
        raise ValueError("Base and noise summaries contain different tools")

    original_non_noise = base.drop(columns=NOISE_COLUMNS).copy(deep=True)
    base.loc[:, NOISE_COLUMNS] = noise.reindex(base.index)[NOISE_COLUMNS]
    pd.testing.assert_frame_equal(
        base.drop(columns=NOISE_COLUMNS),
        original_non_noise,
        check_exact=True,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    base.reset_index().to_csv(args.output, index=False)
    print(f"Updated only {NOISE_COLUMNS} in {args.output}")


if __name__ == "__main__":
    main()
