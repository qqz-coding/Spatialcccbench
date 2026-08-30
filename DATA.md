# Data policy

The GitHub-ready project does not duplicate the full local workspace. Full tool result trees, caches, and intermediate pickles remain external. Each Figure shortcut accepts `--data-root` for those inputs.

The stLearn DLPFC h5ad files used to audit Figure3/Figure5 boundary adaptation are the exception. The original, edge-equalized, gradient-changed, and four signal-loss levels are stored under `result/` with Git LFS. Their checksums are recorded in `result/stlearn_h5ad_manifest.csv`.

Small precomputed CDF tables and the current processed noise/result-summary tables are included where they are directly needed by the published Figure examples.

The spatial-autocorrelation redraws read existing local Moran/Geary arrays from h5ad/CSV files; they do not recompute spatial autocorrelation or permutation statistics.

Upstream tool reproduction is explicit and separate from Figure redraws. `tools/run.py` accepts a local h5ad and CellChat-style LR CSV, writes the established `result/<dataset>/<tool>/` contract, and does not fetch benchmark datasets. SpaTalk and Giotto use spot-by-gene `counts.csv` plus `metadata.csv`, generated from h5ad by `tools/export_h5ad_for_r.py`.
