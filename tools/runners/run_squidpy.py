"""Run the Squidpy ligand-receptor permutation workflow used in the benchmark."""

from __future__ import annotations

import argparse
import time

import numpy as np
import scanpy as sc
import squidpy as sq

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
    parser.add_argument("--n_perms", type=int, default=1000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    np.random.seed(args.random_seed)

    adata = sc.read_h5ad(args.adata_path)
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    adata = filter_housekeeping_genes(adata)

    if args.cluster_key not in adata.obs:
        raise KeyError(f"adata.obs does not contain cluster key '{args.cluster_key}'.")

    lr_frame = load_lr_reference(args.LR_ref_path)
    lr_pairs = list(
        dict.fromkeys(
            (str(ligand), str(receptor))
            for ligand, receptor in zip(lr_frame["ligand.symbol"], lr_frame["receptor.symbol"])
            if str(ligand) in adata.var_names and str(receptor) in adata.var_names
        )
    )
    if not lr_pairs:
        raise ValueError("No ligand-receptor pairs overlap the input genes.")

    started = time.time()
    sq.gr.ligrec(
        adata,
        n_perms=args.n_perms,
        cluster_key=args.cluster_key,
        use_raw=False,
        interactions=lr_pairs,
        seed=args.random_seed,
    )
    result = adata.uns[f"{args.cluster_key}_ligrec"]["pvalues"]
    output_path = output_directory(args, "cellphoneDB") / "result.csv"
    result.to_csv(output_path)
    report_resources(time.time() - started)
    print(f"output: {output_path}")


if __name__ == "__main__":
    main()
