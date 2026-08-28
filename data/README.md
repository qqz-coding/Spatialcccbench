# Optional local data

Only small precomputed tables required by the Figure2 CDF example are included here. Large h5ad objects and full tool results are intentionally external because GitHub has per-file and repository-size limits.

Use `--data-root <path>` with a local data directory that follows the original `dataset/`, `result/`, and precomputed-table layout. Spatial-autocorrelation figures read existing arrays from those files and do not recalculate Moran, Geary, or permutation statistics.
