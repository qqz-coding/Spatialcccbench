#!/usr/bin/env python3
"""Export an h5ad to the spot-by-gene count and metadata CSV contract used by R runners."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--cluster-key", default="cell_type")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    if args.cluster_key not in adata.obs:
        raise KeyError(f"adata.obs does not contain '{args.cluster_key}'.")
    if "spatial" not in adata.obsm:
        raise KeyError("adata.obsm does not contain 'spatial'.")

    matrix = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    counts = pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
    counts.index.name = "spot"
    count_path = output_dir / "counts.csv"
    counts.to_csv(count_path)

    coordinates = np.asarray(adata.obsm["spatial"])
    metadata = pd.DataFrame(
        {
            "spatial1": coordinates[:, 0],
            "spatial2": coordinates[:, 1],
            "cell_type": adata.obs[args.cluster_key].astype(str).to_numpy(),
        },
        index=adata.obs_names,
    )
    metadata.index.name = "spot"
    metadata_path = output_dir / "metadata.csv"
    metadata.to_csv(metadata_path)
    print(f"counts: {count_path}")
    print(f"metadata: {metadata_path}")


if __name__ == "__main__":
    main()
