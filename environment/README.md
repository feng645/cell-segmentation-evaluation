# Software environments

Exports of the environments used in the study. Each segmentation method ran in its own conda
environment; the evaluation pipeline ran in R.

| File | Used for |
|---|---|
| `cellpose.yml` | Cellpose 4.0.6 (Cellpose-SAM); also the ROI crop and Cellpose post-processing scripts |
| `bidcell.yml` | BIDCell 1.0.3, the BIDCell reference scripts, and transcript cleaning |
| `segger.yml` | Segger (see below) |
| `proseg.yml` | Python tools around Proseg (Proseg itself is a Rust program, see below) |
| `baysor_julia/` | `Project.toml` and `Manifest.toml` of the Baysor 0.7.1 Julia environment (Julia 1.10.10) |
| `R_packages.csv` | R and the R packages used by the evaluation pipeline, with versions and sources |
| `segger_dev_b2d2095_local_changes.patch` | Local changes to `segger_dev` that were present when Segger was run |

The conda files were exported with `conda env export --no-builds` on Linux (x86_64), so some
pins may need relaxing on other platforms.

## Python environments

```bash
conda env create -f environment/cellpose.yml
conda env create -f environment/bidcell.yml
conda env create -f environment/proseg.yml
conda env create -f environment/segger.yml
```

## Segger

Segger is not installed from PyPI. It was run from a clone of `segger_dev` at commit `b2d2095`
with local changes (compatibility fixes for the MERSCOPE input and for tiles without training
edges, and model initialization in `train_model.py`):

```bash
git clone https://github.com/EliHei2/segger_dev.git
cd segger_dev
git checkout b2d2095
git apply ../cell-segmentation-evaluation/environment/segger_dev_b2d2095_local_changes.patch
conda activate segger
pip install -e .
```

## Proseg

Proseg 3.0.6 was installed from crates.io:

```bash
cargo install proseg --version 3.0.6
```

## Baysor

Baysor 0.7.1 was installed into Julia 1.10.10 and run through its command-line wrapper with a
precompiled system image. `baysor_julia/Manifest.toml` records the versions of all Julia
dependencies; the `path` entry for Baysor points to the installed package copy on the original
machine. To reinstall:

```julia
using Pkg
Pkg.add(name = "Baysor", version = "0.7.1")
Pkg.build("Baysor")
```

The system image is CPU-specific; it must be built on the same CPU type as the compute nodes.

## R

R 4.4.0 with the packages in `R_packages.csv`. CRAN packages can be installed with
`install.packages()`, Bioconductor packages (release 3.20) with `BiocManager::install()`, and
Cell Marker Accordion from GitHub:

```r
remotes::install_github("TebaldiLab/cellmarkeraccordion", ref = "5ace35d")
```
