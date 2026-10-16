#!/bin/bash
#SBATCH --job-name=segger_boundaries_xenium
#SBATCH -p msismall
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --nice=50
#SBATCH --chdir=/users/9/feng0421/Code/Segger

set -euo pipefail

TRANSCRIPTS_PARQUET="/projects/standard/zhangl4/shared/Data_feng0421/Xenium_segger/predictions/lung_cancer_v1_pred_0.5_True_4_6.0_5_5.0_20260406/segger_transcripts.parquet"
GEOJSON_OUTPUT="/projects/standard/zhangl4/shared/Data_feng0421/Xenium_segger/predictions/lung_cancer_v1_pred_0.5_True_4_6.0_5_5.0_20260406/segger_cell_boundaries.geojson.gz"

echo "[info] job_id=${SLURM_JOB_ID}"
echo "[info] hostname=$(hostname)"
echo "[info] start_time=$(date -Iseconds)"
echo "[info] transcripts=${TRANSCRIPTS_PARQUET}"
echo "[info] output=${GEOJSON_OUTPUT}"

module purge
module load miniforge gcc clang
eval "$(conda shell.bash hook)"
conda activate segger

python /users/9/feng0421/Code/Segger/segger_transcripts_to_boundaries.py \
  --transcripts-parquet "${TRANSCRIPTS_PARQUET}" \
  --geojson-output      "${GEOJSON_OUTPUT}" \
  --cell-id-col         segger_cell_id \
  --min-points          3

echo "[info] end_time=$(date -Iseconds)"
echo "[info] output written to: ${GEOJSON_OUTPUT}"
