# Tutorials

The root `tutorial.ipynb` is the concise entry point. Reusable tutorial code is kept under `tutorials/scripts/`:

- `benchmark_accuracy.py`: extract existing CCC results and calculate precision/recall/F1.
- `redraw_precomputed_spatial_figures.py`: redraw Figure2 or FigureS3 from existing h5ad/CSV arrays.
- `redraw_noise_robustness.py`: redraw Figure4 from processed noise metrics.
- `redraw_boundary_adaptation.py`: redraw Figure3 from existing perturbation results.

The Figure and evaluation scripts require local data paths and never launch upstream CCC tools implicitly. Explicit environment installation and upstream execution are handled separately by `tools/run.py`; see `tools/README.md` and `tools/tutorial.ipynb`.
