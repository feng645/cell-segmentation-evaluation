#!/bin/bash
#SBATCH --job-name=cellpose_pipeline
#SBATCH --cpus-per-task=8
#SBATCH -p msigpu
#SBATCH --gres=gpu:1
#SBATCH --time=2:00:00
#SBATCH --mem=64GB
#SBATCH --chdir=/users/9/feng0421/Code/Cellpose

set -euo pipefail

DAPI_PATH="/users/9/feng0421/BIDCell/data/MerFish/roi_dapi.tif"
TRANSCRIPTS_PATH="/users/9/feng0421/BIDCell/data/MerFish/roi_transcripts.csv"
AFFINE_MATRIX="/users/9/feng0421/BIDCell/data/MerFish/roi_micron_to_pixel_transform.csv"
ROI_UM=""
OUTPUT_DIR="/users/9/feng0421/Code/R/Evaluation/Cellpose_out"
EXPAND_UM="5"
PLATFORM="merfish"
RUN_PLOT="true"

usage() {
  cat <<'EOF'
Usage:
  sbatch run_cellpose_pipeline.sh [--roi-um X_MIN Y_MIN X_MAX Y_MAX | --affine-matrix PATH]

Optional arguments:
  --dapi PATH
  --transcripts PATH
  --output-dir PATH
  --expand-um FLOAT
  --platform {xenium|cosmx|merfish}
  --no-plot

Defaults:
  --dapi /users/9/feng0421/BIDCell/data/MerFish/roi_dapi.tif
  --transcripts /users/9/feng0421/BIDCell/data/MerFish/roi_transcripts.csv
  --affine-matrix /users/9/feng0421/BIDCell/data/MerFish/roi_micron_to_pixel_transform.csv
EOF
}

COORD_MODE_SET="false"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dapi)
      DAPI_PATH="$2"
      shift 2
      ;;
    --transcripts)
      TRANSCRIPTS_PATH="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --expand-um)
      EXPAND_UM="$2"
      shift 2
      ;;
    --platform)
      PLATFORM="$2"
      shift 2
      ;;
    --affine-matrix)
      if [[ "$COORD_MODE_SET" == "true" ]]; then
        echo "Only one of --roi-um or --affine-matrix may be provided." >&2
        exit 1
      fi
      AFFINE_MATRIX="$2"
      ROI_UM=""
      COORD_MODE_SET="true"
      shift 2
      ;;
    --roi-um)
      if [[ "$COORD_MODE_SET" == "true" ]]; then
        echo "Only one of --roi-um or --affine-matrix may be provided." >&2
        exit 1
      fi
      ROI_UM="$2 $3 $4 $5"
      AFFINE_MATRIX=""
      COORD_MODE_SET="true"
      shift 5
      ;;
    --no-plot)
      RUN_PLOT="false"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

mkdir -p "$OUTPUT_DIR"

MASK_GZ="$OUTPUT_DIR/cellpose-masks_roi.npy.gz"
CELLPROB_GZ="$OUTPUT_DIR/cellpose-cellprobs_roi.npy.gz"
EXPANDED_MASK_GZ="$OUTPUT_DIR/cellpose-masks_roi_expanded_${EXPAND_UM}um.npy.gz"
EXPANDED_TIF="$OUTPUT_DIR/cellpose_masks_roi_expanded_${EXPAND_UM}um.tif"
CELLPOSE_RESOURCE_LOG="$OUTPUT_DIR/cellpose_step_resources.log"
CELLPOSE_PEAK_MEM_FILE="$OUTPUT_DIR/cellpose_step_peak_mem.txt"

COORD_ARGS=()
if [[ -n "$ROI_UM" ]]; then
  read -r XMIN YMIN XMAX YMAX <<< "$ROI_UM"
  COORD_ARGS=(--roi-um "$XMIN" "$YMIN" "$XMAX" "$YMAX")
else
  COORD_ARGS=(--affine-matrix "$AFFINE_MATRIX")
fi

echo "Cellpose pipeline started: $(date)"

module load miniforge
eval "$(conda shell.bash hook)"
conda activate cellpose

CELLPOSE_ARGS=(
  python Cellpose.py
  "$DAPI_PATH"
  --platform "$PLATFORM"
  --masks-output "$MASK_GZ"
  --cellprobs-output "$CELLPROB_GZ"
)
if [[ "$RUN_PLOT" == "true" ]]; then
  CELLPOSE_ARGS+=(--plot)
fi

echo "Writing Cellpose step resource log to: $CELLPOSE_RESOURCE_LOG"
/usr/bin/time -v "${CELLPOSE_ARGS[@]}" 2> "$CELLPOSE_RESOURCE_LOG"
grep "Maximum resident set size" "$CELLPOSE_RESOURCE_LOG" > "$CELLPOSE_PEAK_MEM_FILE"
echo "Writing Cellpose peak memory summary to: $CELLPOSE_PEAK_MEM_FILE"

python expand_cellpose_labels.py \
  --mask-gz "$MASK_GZ" \
  "${COORD_ARGS[@]}" \
  --expand-um "$EXPAND_UM" \
  --out-mask-gz "$EXPANDED_MASK_GZ" \
  --out-label-tif "$EXPANDED_TIF"

python prepare_cellpose_for_evaluation.py \
  --mask-gz "$EXPANDED_MASK_GZ" \
  --transcripts "$TRANSCRIPTS_PATH" \
  "${COORD_ARGS[@]}" \
  --out-dir "$OUTPUT_DIR" \
  --out-label-tif "$EXPANDED_TIF"

module load R/4.4.0-openblas-rocky8

cd /users/9/feng0421/Code/R/Evaluation/Main/Morphology
Rscript run_cellpose_morphology.R

cd /users/9/feng0421/Code/R/Evaluation/Main/Expression
Rscript assign_transcripts_to_cells.R
Rscript run_expression_eval.R

cd /users/9/feng0421/Code/R/Evaluation/Main/Spatial
Rscript run_spatial_eval.R

cd /users/9/feng0421/Code/R/Evaluation/Main
Rscript BANKSY.R

echo "Cellpose pipeline finished: $(date)"
