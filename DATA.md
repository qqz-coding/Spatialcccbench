# Data policy

The GitHub-ready project does not duplicate the full local workspace. Large h5ad files, full tool result trees, caches, and intermediate pickles remain external. Each Figure shortcut accepts `--data-root` for those inputs.

Small precomputed CDF tables and the current processed noise/result-summary tables are included where they are directly needed by the published Figure examples.

The spatial-autocorrelation redraws read existing local Moran/Geary arrays from h5ad/CSV files; they do not recompute spatial autocorrelation or permutation statistics.
