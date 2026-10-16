# MainV2 Scripts

This directory contains the runnable entry points for the `MainV2` evaluation pipeline.

## Recommended execution order

1. `run_phase1_eval.R`
   Builds the unified transcript universe, standardizes method boundaries, assigns transcripts to cells, and writes method-level count tables.
2. `run_annotation_eval.R`
   Runs Seurat-based QC, clustering, and marker-driven annotation on unified per-method count tables.
3. `run_expression_eval.R`
   Computes marker purity, contamination, and count summaries using the `annotation/` outputs.
4. `run_morphology_eval.R`
   Computes polygon-derived geometric metrics and summarizes them by annotated cell type.
5. `run_spatial_eval.R`
   Computes KNN-based spatial mixing, domain-backend internal consistency scores (default `BANKSY`, optional `script` backend), and neighborhood marker leakage using morphology centroids plus assignment-derived expression counts.
6. `run_visualization_eval.R`
   Produces method-level visualization panels for assignments, QC, centroids, boundaries, and zoom-ins.

## Domain Backends

- Default: `spatial$domain_clustering$backend = "banksy"`
  Runs the generic BANKSY backend on standardized per-method inputs (`cells_eval.csv` + assignment-derived counts).
- Optional: `spatial$domain_clustering$backend = "script"`
  Calls a user-provided script with positional arguments
  `<input_manifest.json> <output_labels.csv>`.
  The script must write a CSV containing at least `cell_id` and `domain_label`.

Example wrapper:

```bash
bash /users/9/feng0421/Code/R/Evaluation/MainV2/Scripts/example_domain_backend_banksy.sh \
  /path/to/input_manifest.json \
  /path/to/output_labels.csv
```


## One-shot runner

- `run_all_eval.R`
  Runs phase1, annotation, expression, morphology, spatial, and visualization in sequence.


## Cluster submission

- `submit_xenium_lung_ffpe_mainv2.sh`
  Example SLURM submission script for the Xenium run configuration that writes outputs to the shared project space.

## Usage

All R entry points accept either a run config or, where noted, an `output_root`.

Example:

```bash
module load R/4.4.0-openblas-rocky8
Rscript /users/9/feng0421/Code/R/Evaluation/MainV2/Scripts/run_all_eval.R \
  /users/9/feng0421/Code/R/Evaluation/MainV2/Config/Runs/merscope_lung_roi.R
```
