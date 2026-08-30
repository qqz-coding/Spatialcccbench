# Environment definitions

These are minimal, named Conda definitions derived from the manuscript machine snapshots. They intentionally omit Linux build hashes, absolute `prefix:` paths, university mirror URLs, Jupyter front-end packages, and unrelated transitive packages from the original full exports.

Each YAML includes `nodefaults`, preventing a user's global `.condarc` mirrors or default channels from silently changing the solve.

| Environment | Historical core versions |
|---|---|
| `spatialccc-squidpy` | Python 3.9, AnnData 0.10.9, Scanpy 1.10.0, Squidpy 1.2.2 |
| `spatialccc-commot` | Python 3.11, AnnData 0.11.3, Scanpy 1.11.0, COMMOT 0.0.3 |
| `spatialccc-spatialdm` | Python 3.9, AnnData 0.10.x, Scanpy 1.10.3, SpatialDM 0.2.0 |
| `spatialccc-stlearn` | Python 3.8, stLearn 0.4.12, TensorFlow 2.4.1, Scanpy 1.9.8 |
| `spatialccc-spatalk` | R 4.2 plus SpaTalk and NNLM from their upstream repositories |
| `spatialccc-giotto` | R 4.2 plus Giotto v3.3.2 |

Install through the dispatcher rather than calling the YAML files manually:

```bash
python tools/run.py install squidpy commot spatialdm
python tools/run.py install stlearn
python tools/run.py install spatalk giotto
```

`install_r_packages.R` is called automatically after the R Conda environments are created. Network access and a working compiler toolchain may be required. stLearn 0.4.12 and TensorFlow 2.4.1 are most reliably reproduced on Linux or WSL.

The environment files reproduce the recorded direct versions but are not byte-for-byte lock files for every operating system. After installation, keep `conda list --explicit` and `sessionInfo()` with the run log when archival-level provenance is required.
