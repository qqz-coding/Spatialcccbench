# Reproducing benchmark tool runs

`tools/run.py` is the single entry point for environment setup and tool execution. Each upstream tool remains in an isolated Conda environment because the manuscript versions have incompatible Python, NumPy, TensorFlow, and R requirements.

## Quick start

Run these commands from the repository root:

```bash
python tools/run.py list
python tools/run.py doctor
python tools/run.py install squidpy
python tools/run.py run squidpy \
  --dataset st_lymphnode \
  --adata data/st_lymphnode.h5ad \
  --lr-db data/cellchatDB_human.csv \
  --output-root result
```

The dispatcher uses `micromamba`, `mamba`, or `conda`, in that order. Override discovery with `--manager conda`. It calls `conda run`, so shell activation is not required.

For a one-command setup and run, replace the separate `install` and `run` calls with:

```bash
python tools/run.py reproduce squidpy \
  --dataset st_lymphnode \
  --adata data/st_lymphnode.h5ad \
  --lr-db data/cellchatDB_human.csv
```

Preview any install or run without changing the machine:

```bash
python tools/run.py install squidpy commot spatialdm --dry-run
python tools/run.py run commot --dataset demo --adata demo.h5ad \
  --lr-db cellchatDB_human.csv --dry-run
```

## Tools and outputs

| Tool ID | Environment | Output below `result/<dataset>/` | State |
|---|---|---|---|
| `squidpy` | `spatialccc-squidpy` | `cellphoneDB/result.csv` | Ready |
| `baseline_1` | `spatialccc-squidpy` | `baseline_1/result.csv` | Ready |
| `baseline_2` | `spatialccc-squidpy` | `baseline_2/result.csv` | Ready |
| `commot` | `spatialccc-commot` | `COMMOT/result.h5ad` | Ready |
| `spatialdm` | `spatialccc-spatialdm` | `spatialDM/result.csv` | Ready |
| `stlearn` | `spatialccc-stlearn` | `stlearn/result.h5ad` | Ready, Linux/WSL recommended |
| `stlearn_without_spotmixture` | `spatialccc-stlearn` | `stlearn/result_without_mixture.h5ad` | Ready, Linux/WSL recommended |
| `spatalk` | `spatialccc-spatalk` | `spatalk/result.csv` | Manuscript-era R workflow |
| `giotto` | `spatialccc-giotto` | `giotto/result.csv` | Giotto `v3.3.2` workflow |
| `cellagentchat` | External | `cellagentchat/result.csv` | Upstream runner unavailable |

The Squidpy result directory remains named `cellphoneDB` to preserve the published `toolkit.preprocess` contract. This does not invoke the CellPhoneDB package.

## Inputs

Python runners use:

- `--adata`: h5ad with expression in `X`, coordinates in `obsm['spatial']`, and cell types in `obs['cell_type']`.
- `--lr-db`: CSV containing `ligand.symbol` and `receptor.symbol`. `pathway_name` and `annotation` are also used where required.
- `--dataset`: output label, including nested labels such as `noise/dropout`.

Pass a runner-specific option as one token per `--runner-arg`:

```bash
python tools/run.py run squidpy --dataset demo --adata demo.h5ad \
  --lr-db lr.csv --runner-arg=--n_perms --runner-arg=100
```

For SpaTalk and Giotto, first export the common R CSV contract:

```bash
conda run -n spatialccc-squidpy python tools/export_h5ad_for_r.py \
  --adata data/st_lymphnode.h5ad \
  --output-dir data/r_export/st_lymphnode

python tools/run.py run spatalk \
  --dataset st_lymphnode \
  --count-csv data/r_export/st_lymphnode/counts.csv \
  --meta-csv data/r_export/st_lymphnode/metadata.csv \
  --lr-db data/cellchatDB_human.csv
```

`counts.csv` is spot by gene. `metadata.csv` has matching spot IDs plus `spatial1`, `spatial2`, and `cell_type`.

## Batch execution

Copy `tools/jobs.example.csv`, edit paths, then run:

```bash
python tools/run.py batch --jobs tools/jobs.example.csv
```

Jobs run sequentially. Each run writes a combined stdout/stderr log to `result/<dataset>/_logs/<tool>.log` and checks that the expected result exists before reporting success.

## Reproducibility boundary

The runners preserve manuscript preprocessing, permutation counts, output schemas, and stLearn spot-mixture modes. They execute CCC tools but do not launch the separate post hoc Moran/Geary workflow used for the spatial-autocorrelation Figures.

Every archived `run_CellAgentChat.py` found in the deposited project is an empty 0-byte file. The dispatcher therefore reports CellAgentChat as unavailable rather than inventing an implementation. Existing CellAgentChat result CSV files remain readable by `toolkit.preprocess`.

Upstream installation references: [Squidpy](https://squidpy.readthedocs.io/), [COMMOT 0.0.3](https://commot.readthedocs.io/en/latest/installation.html), [SpatialDM](https://spatialdm.readthedocs.io/en/latest/install.html), [stLearn](https://stlearn.readthedocs.io/), [SpaTalk](https://github.com/ZJUFanLab/SpaTalk), and [Giotto tags](https://github.com/giotto-suite/Giotto/tags).
