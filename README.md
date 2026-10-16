# cell-segmentation-evaluation

Code, run configurations, platform profiles, and scripts for the study *Cell Segmentation Choice
Shapes Downstream Analyses in Imaging-Based Spatial Transcriptomics*. The study compares four
hybrid cell segmentation methods (BIDCell, Baysor, Proseg, Segger) and an image-only baseline
(Cellpose) on a MERSCOPE lung cancer region of interest, and applies the same evaluation workflow
to a 10x Xenium lung FFPE section.

The repository has two parts:

1. **Segmentation runs** (`Scripts/segmentation/`): the commands, configuration files, and Slurm
   scripts used to run each method.
2. **Evaluation workflow** (`Config/`, `Source/`, `Scripts/`, `Extra/`): a platform-agnostic R
   pipeline that reads every method's cell boundaries against one canonical transcript table,
   reassigns transcripts, and computes morphology, expression, and spatial-context metrics with
   the same QC and annotation for all methods.

## Layout

```
Config/
  Profiles/          platform profiles (column conventions, QC thresholds, defaults)
  Runs/              run configurations used in the study
Source/              evaluation pipeline (R; run_phase1_eval.py is an optional Python phase 1)
Scripts/             pipeline entry points and Slurm submission scripts (see Scripts/README.md)
  segmentation/
    data_preparation/  MERSCOPE ROI crop and transcript cleaning
    cellpose/          nuclear segmentation (Cellpose-SAM) and 5 um expansion
    bidcell/           BIDCell configs and run script; reference/ builds the single-cell reference
    baysor/            Baysor configs and Slurm scripts
    proseg/            Proseg run script
    segger/            Segger data preparation, training, and inference scripts
Extra/               cell-type marker tables used for annotation and marker purity
```

## Data

| Data | Source |
|---|---|
| MERSCOPE FFPE human lung cancer | Vizgen Data Release Program, https://info.vizgen.com/ffpe-showcase |
| Xenium FFPE human lung cancer | 10x Genomics, https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard |
| Human Lung Cell Atlas (BIDCell and Segger reference) | CELLxGENE, https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293 |
| NSCLC single-cell data, GSE148071 (BIDCell reference) | GEO accession GSE148071 |

The MERSCOPE analysis uses the region x, y in [550, 800] um (global coordinates).

## Workflow

1. **Prepare the transcripts** (`Scripts/segmentation/data_preparation/`).
   `run_roi_merfish.sh` crops the MERSCOPE ROI (transcripts and DAPI) with `roi_crop_merfish.py`.
   `clean_trans.py` removes blank and control codewords and the genes listed in
   `gene_to_filter.txt` (MERSCOPE), giving the canonical transcript table;
   `run_clean_xenium_transcripts.sh` does the same for Xenium.
2. **Run the segmentation methods** (`Scripts/segmentation/<method>/`). Each folder holds the
   script that was submitted. Inputs, priors, and non-default settings are summarized in
   Online Resource 1, Table S4 of the paper.
3. **Run the evaluation**:

   ```bash
   module load R/4.4.0-openblas-rocky8 gdal/3.8.2-rocky8 geos/3.11.2-rocky8 udunits/2.2.27.6_gcc11.3.0-rocky8
   Rscript Scripts/run_all_eval.R Config/Runs/merscope_lung_roi_gene_rate_baysorfix.R
   ```

   `Scripts/submit_merscope_lung_roi_baysorfix.sh` and
   `Scripts/submit_xenium_lung_ffpe_baysorfix.sh` are the Slurm scripts for the reported results.
   Stages can also be run one at a time (see `Scripts/README.md`).

## Run configurations

| Configuration | Use |
|---|---|
| `merscope_lung_roi_gene_rate_baysorfix.R`, `xenium_lung_ffpe_shared_gene_rate_baysorfix.R` | Reported results: abundance-normalized marker scaling, Baysor rerun with the correct unassigned prior label |
| `merscope_lung_roi_gene_rate.R`, `xenium_lung_ffpe_shared_gene_rate.R` | Same, with the earlier Baysor run |
| `merscope_lung_roi.R`, `xenium_lung_ffpe_shared.R` | Base configurations the others build on |
| `*_banksy_backend.R`, `*_script_backend.R` | Examples of the two spatial-domain backends |

## Notes

- **Paths.** Configurations and scripts contain the absolute paths of the cluster where the study
  was run. Edit `canonical_transcripts`, the method `input$path` entries, `markers_path`, and
  `output_root` in the run configuration before reusing it.
- **Baysor.** The reported Baysor results come from `Scripts/segmentation/baysor/`, which sets
  `unassigned_prior_label` to the vendor value (`-1` for MERSCOPE, `UNASSIGNED` for Xenium). Baysor
  0.7.1 has no command-line flag for this option, so it is set in the TOML config.
- **Cellpose.** `Cellpose.py` was run with Cellpose 4.0.6, which uses the Cellpose-SAM model and
  ignores the `model_type` and `channels` arguments.
- **Proseg.** `run_proseg.sh` contains the Xenium command; the MERSCOPE command (`--merscope`,
  `--voxel-layers 7`) is kept below it as a comment.
- **Segger** was run from a clone of `segger_dev` (version 0.1.0, commit `b2d2095`).

## Software

| Software | Version |
|---|---|
| Cellpose | 4.0.6 |
| BIDCell | 1.0.3 |
| Baysor | 0.7.1 (Julia 1.10.10) |
| Proseg | 3.0.6 |
| Segger | 0.1.0 (`segger_dev` commit `b2d2095`; PyTorch 2.5.0, PyG 2.7.0) |
| R | 4.4.0 |
| Seurat / SeuratObject | 5.4.0 / 5.3.0 |
| Cell Marker Accordion | 1.0.0 |
| Banksy | 1.2.0 |
| SpatialExperiment | 1.16.0 |
| FNN | 1.1.4.1 |
| uwot | 0.2.4 |
