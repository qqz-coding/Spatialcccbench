"""Shared stLearn workflow for the two spot-mixture modes."""

from __future__ import annotations

import argparse
import time

import numpy as np
import scanpy as sc
import stlearn as st

from _common import (
    add_common_arguments,
    filter_housekeeping_genes,
    get_cell_type_matrix,
    load_lr_reference,
    output_directory,
    report_resources,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    parser.add_argument("--n_pairs", type=int, default=1000)
    parser.add_argument("--n_perms", type=int, default=1000)
    parser.add_argument("--n_cpus", type=int, default=2)
    parser.add_argument("--min_spots", type=int, default=3)
    parser.add_argument("--cell_prop_cutoff", type=float, default=0.2)
    parser.add_argument("--pval_adj_cutoff", type=float, default=0.05)
    parser.add_argument("--spatial_library_key", default="auto")
    return parser.parse_args()


def prepare_lr_pairs(path: str) -> np.ndarray:
    frame = load_lr_reference(path)
    pairs: list[str] = []
    for ligand, receptors in zip(frame["ligand.symbol"], frame["receptor.symbol"]):
        for receptor in str(receptors).split(", "):
            pairs.append(f"{ligand}_{receptor}")
    return np.unique(np.asarray(pairs, dtype=str))


def set_spatial_quality(adata, requested_key: str) -> None:
    spatial = adata.uns.get("spatial")
    if not isinstance(spatial, dict) or not spatial:
        return
    key = next(iter(spatial)) if requested_key == "auto" else requested_key
    if key not in spatial:
        raise KeyError(f"Spatial library '{key}' is absent from adata.uns['spatial'].")
    if isinstance(spatial[key], dict):
        spatial[key]["use_quality"] = "hires"


def main(*, spot_mixtures: bool) -> None:
    args = parse_args()
    np.random.seed(args.random_seed)
    lr_pairs = prepare_lr_pairs(args.LR_ref_path)

    adata = sc.read_h5ad(args.adata_path)
    adata.var_names_make_unique()
    set_spatial_quality(adata, args.spatial_library_key)
    if spot_mixtures:
        adata.uns[args.cluster_key] = get_cell_type_matrix(adata, args.cluster_key)
    elif args.cluster_key not in adata.obs:
        raise KeyError(f"adata.obs does not contain cluster key '{args.cluster_key}'.")

    st.pp.filter_genes(adata, min_cells=1)
    st.pp.normalize_total(adata)
    adata = filter_housekeeping_genes(adata)
    adata.var_names = adata.var_names.astype(str).str.replace("_", "-", regex=False)
    result_adata = adata.copy()

    started = time.time()
    st.tl.cci.run(
        adata,
        lr_pairs,
        min_spots=args.min_spots,
        distance=None,
        n_pairs=args.n_pairs,
        n_cpus=args.n_cpus,
    )
    st.tl.cci.adj_pvals(
        adata,
        correct_axis="spot",
        pval_adj_cutoff=args.pval_adj_cutoff,
        adj_method="fdr_bh",
    )
    st.tl.cci.run_cci(
        adata,
        args.cluster_key,
        min_spots=args.min_spots,
        spot_mixtures=spot_mixtures,
        n_cpus=args.n_cpus,
        cell_prop_cutoff=args.cell_prop_cutoff,
        sig_spots=True,
        n_perms=args.n_perms,
    )
    result_adata.uns["per_lr_cci_pvals_cell_type"] = adata.uns["per_lr_cci_pvals_cell_type"]
    filename = "result.h5ad" if spot_mixtures else "result_without_mixture.h5ad"
    output_path = output_directory(args, "stlearn") / filename
    result_adata.write_h5ad(output_path)
    report_resources(time.time() - started)
    print(f"output: {output_path}")
