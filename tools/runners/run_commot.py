"""Run the COMMOT cluster-level spatial permutation workflow."""

from __future__ import annotations

import argparse
import time

import commot as ct
import pandas as pd
import scanpy as sc

from _common import (
    add_common_arguments,
    filter_housekeeping_genes,
    load_lr_reference,
    output_directory,
    report_resources,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    parser.add_argument("--database_name", default="lymnode")
    parser.add_argument("--distance_threshold", type=float, default=500.0)
    parser.add_argument("--min_cell_pct", type=float, default=0.05)
    return parser.parse_args()


def commot_database(path: str) -> pd.DataFrame:
    frame = load_lr_reference(path).copy()
    pathway = frame["pathway_name"] if "pathway_name" in frame else "custom"
    receptors = frame["receptor.symbol"].astype(str).str.replace(r",\s*", "_", regex=True)
    return pd.DataFrame(
        {
            0: frame["ligand.symbol"].astype(str),
            1: receptors,
            2: pathway,
        }
    ).drop_duplicates()


def main() -> None:
    args = parse_args()
    adata = sc.read_h5ad(args.adata_path)
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    adata = filter_housekeeping_genes(adata)
    if args.cluster_key not in adata.obs:
        raise KeyError(f"adata.obs does not contain cluster key '{args.cluster_key}'.")

    database = commot_database(args.LR_ref_path)
    database = ct.pp.filter_lr_database(database, adata, min_cell_pct=args.min_cell_pct)
    if database.empty:
        raise ValueError("No ligand-receptor pairs passed COMMOT filtering.")

    started = time.time()
    ct.tl.cluster_communication_spatial_permutation(
        adata,
        database,
        database_name=args.database_name,
        clustering=args.cluster_key,
        dis_thr=args.distance_threshold,
    )
    output_path = output_directory(args, "COMMOT") / "result.h5ad"
    adata.write_h5ad(output_path)
    report_resources(time.time() - started)
    print(f"output: {output_path}")


if __name__ == "__main__":
    main()
