# Optional local data

Only small precomputed tables required by the Figure2 CDF example are included here. Full tool result trees remain external. The stLearn DLPFC h5ad inputs used for the Figure3/Figure5 boundary audit are versioned separately under `../result/` with Git LFS.

Use `--data-root <path>` with a local data directory that follows the original `dataset/`, `result/`, and precomputed-table layout. Spatial-autocorrelation figures read existing arrays from those files and do not recalculate Moran, Geary, or permutation statistics.
