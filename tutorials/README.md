# Tutorials

The root `tutorial.ipynb` is the concise entry point. Reusable tutorial code is kept under `tutorials/scripts/`:

- `benchmark_accuracy.py`: extract existing CCC results and calculate precision/recall/F1.
- `redraw_precomputed_spatial_figures.py`: redraw Figure2 or FigureS3 from existing h5ad/CSV arrays.
- `redraw_noise_robustness.py`: redraw Figure4 from processed noise metrics.
- `redraw_boundary_adaptation.py`: redraw Figure3 from existing perturbation results.

The scripts require local data paths and never launch the upstream CCC tools implicitly.
