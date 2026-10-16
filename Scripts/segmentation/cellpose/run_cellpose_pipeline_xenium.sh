#!/bin/bash
#SBATCH --job-name=cellpose_xenium
#SBATCH --cpus-per-task=8
#SBATCH -p msigpu
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --chdir=/users/9/feng0421/Code/Cellpose
set -euo pipefail

DATA_DIR="/projects/standard/zhangl4/shared/Data_feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs"
DAPI_PATH="$DATA_DIR/morphology_focus_dapi_2d.tif"
OME_PATH="$DATA_DIR/morphology.ome.tif"
TRANSCRIPTS_PATH="$DATA_DIR/transcripts_filtered.csv"
OUT_DIR="/projects/standard/zhangl4/shared/Data_feng0421/Xenium_Cellpose"

AFFINE_MATRIX="$OUT_DIR/xenium_micron_to_pixel_transform.csv"
MASK_GZ="$OUT_DIR/cellpose-masks_roi.npy.gz"
CELLPROB_GZ="$OUT_DIR/cellpose-cellprobs_roi.npy.gz"
EXPAND_UM="5"
TILE_SIZE="8192"
EXPANDED_MASK_GZ="$OUT_DIR/cellpose-masks_roi_expanded_${EXPAND_UM}um.npy.gz"
EXPANDED_TIF="$OUT_DIR/cellpose_masks_roi_expanded_${EXPAND_UM}um.tif"
CELLPOSE_RESOURCE_LOG="$OUT_DIR/cellpose_segmentation_time.log"
CELLPOSE_PEAK_MEM_FILE="$OUT_DIR/cellpose_step_peak_mem.txt"

mkdir -p "$OUT_DIR"

echo "Xenium Cellpose job started: $(date)"
echo "DATA_DIR=$DATA_DIR"
echo "DAPI_PATH=$DAPI_PATH"
echo "OME_PATH=$OME_PATH"
echo "TRANSCRIPTS_PATH=$TRANSCRIPTS_PATH"
echo "OUT_DIR=$OUT_DIR"

module load miniforge
eval "$(conda shell.bash hook)"
conda activate cellpose

python generate_xenium_affine.py \
  --ome-tif "$OME_PATH" \
  --out "$AFFINE_MATRIX"

echo "Affine matrix written to $AFFINE_MATRIX"

echo "Running Cellpose with tile size $TILE_SIZE"
{
  echo "Segmentation step started: $(date)"
  echo "Host: $(hostname)"
  echo "Working directory: $(pwd)"
  echo "SLURM job ID: ${SLURM_JOB_ID:-N/A}"
  echo "Timing command: /usr/bin/time -v python Cellpose.py ..."
  echo
} > "$CELLPOSE_RESOURCE_LOG"
/usr/bin/time -v python Cellpose.py \
  "$DAPI_PATH" \
  --platform xenium \
  --tile-size "$TILE_SIZE" \
  --masks-output "$MASK_GZ" \
  --cellprobs-output "$CELLPROB_GZ" \
  2>> "$CELLPOSE_RESOURCE_LOG"

{
  echo
  echo "Segmentation step finished: $(date)"
  echo "Segmentation log saved to: $CELLPOSE_RESOURCE_LOG"
} >> "$CELLPOSE_RESOURCE_LOG"

grep "Maximum resident set size" "$CELLPOSE_RESOURCE_LOG" > "$CELLPOSE_PEAK_MEM_FILE" || true
echo "Peak memory summary saved to: $CELLPOSE_PEAK_MEM_FILE"

python expand_cellpose_labels.py \
  --mask-gz "$MASK_GZ" \
  --affine-matrix "$AFFINE_MATRIX" \
  --expand-um "$EXPAND_UM" \
  --out-mask-gz "$EXPANDED_MASK_GZ" \
  --out-label-tif "$EXPANDED_TIF"

python prepare_cellpose_for_evaluation.py \
  --mask-gz "$EXPANDED_MASK_GZ" \
  --transcripts "$TRANSCRIPTS_PATH" \
  --affine-matrix "$AFFINE_MATRIX" \
  --out-dir "$OUT_DIR" \
  --out-label-tif "$EXPANDED_TIF"

echo "Prepared outputs:"
echo "  $AFFINE_MATRIX"
echo "  $MASK_GZ"
echo "  $CELLPROB_GZ"
echo "  $EXPANDED_MASK_GZ"
echo "  $EXPANDED_TIF"
echo "  $OUT_DIR/cellpose_transcript_assignments.csv"
echo "  $OUT_DIR/cellpose_expr_mat.csv"
echo "  $OUT_DIR/cellpose_centroids.csv"
echo "  $OUT_DIR/cellpose_boundaries.geojson.gz"
echo "  $OUT_DIR/cellpose_mask_stats.csv"

echo "Xenium Cellpose job finished: $(date)"
