# SpatialCCCbench

A GitHub-ready organization of the SpatialCCC benchmark toolkit, tutorial entry points, and manuscript Figure redraws.

## Start here

```bash
python scripts/validate_repository.py
```

- [tutorial.ipynb](tutorial.ipynb): concise tutorial with Markdown guidance and executable entry points.
- [figures_index.csv](figures_index.csv): one-to-one Figure, shortcut, notebook, and result mapping.
- [figures/](figures/): Figure1-5 and FigureS1-S9, each with `run.py`, a concise demo notebook, and `results/published/`.
- [results_manifest.csv](results_manifest.csv): SHA-256 inventory of bundled Figure inputs, submitted SVGs, and current script-generated results.
- [tutorials/scripts/](tutorials/scripts/): extracted tutorial workflows.
- [toolkit/](toolkit/): reusable result harmonization, metrics, plotting, and analysis modules.
- [environments/](environments/): tool-specific Conda definitions.
- [DATA.md](DATA.md): external data policy and `--data-root` convention.

## Figure commands

Quickly inspect a published result:

```bash
python figures/Figure4/run.py
```

Redraw a computational Figure when local inputs are available:

```bash
python figures/Figure2/run.py --reproduce --data-root /path/to/local/data --step 2
python figures/FigureS3/run.py --reproduce --data-root /path/to/local/data
```

Spatial-autocorrelation Figure scripts read precomputed h5ad/CSV arrays. They do not recalculate Moran, Geary, or permutation statistics.

## Repository boundary

The canonical, one-to-one Figure organization is under `figures/`; full local result trees and generated raster caches are excluded by `.gitignore`. Compact SVG/CSV/JSON redraw outputs are versioned. The stLearn DLPFC h5ad inputs required to audit Figure3/Figure5 boundary adaptation are explicitly versioned under `result/` with Git LFS.

The submitted SVG is stored as `results/published/<Figure>.svg`. Computational redraws are written to `results/generated/` and remain separate from the submitted final layout, especially where the final figure contains manual assembly. Figure3 ships a compact ligand-count table exported from the latest raw tool results, so its panels can be redrawn without cloning the multi-gigabyte local result tree.
