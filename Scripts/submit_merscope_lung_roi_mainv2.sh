#!/bin/bash
#SBATCH --job-name=mainv2_merscope_eval
#SBATCH -p msismall
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --nice=50
#SBATCH --array=0-1
#SBATCH --chdir=/users/9/feng0421/Code/R/Evaluation/MainV2

set -euo pipefail

PIPELINE_ROOT="/users/9/feng0421/Code/R/Evaluation/MainV2"
VARIANTS=(
  "banksy:/users/9/feng0421/Code/R/Evaluation/MainV2/Config/Runs/merscope_lung_roi_banksy_backend.R:/users/9/feng0421/Code/R/Evaluation/MainV2/Outputs/merscope_lung_roi_banksy_backend"
  "script:/users/9/feng0421/Code/R/Evaluation/MainV2/Config/Runs/merscope_lung_roi_script_backend.R:/users/9/feng0421/Code/R/Evaluation/MainV2/Outputs/merscope_lung_roi_script_backend"
)

TASK_INDEX="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "${TASK_INDEX}" -lt 0 || "${TASK_INDEX}" -ge "${#VARIANTS[@]}" ]]; then
  echo "[error] Invalid SLURM_ARRAY_TASK_ID=${TASK_INDEX}" >&2
  exit 1
fi

IFS=":" read -r VARIANT_NAME RUN_CONFIG OUTPUT_ROOT <<< "${VARIANTS[${TASK_INDEX}]}"
LOG_DIR="${OUTPUT_ROOT}/logs"

mkdir -p "${OUTPUT_ROOT}" "${LOG_DIR}"

exec > >(tee -a "${LOG_DIR}/mainv2_merscope_eval_${VARIANT_NAME}_${SLURM_JOB_ID}_${TASK_INDEX}.log")
exec 2>&1

echo "[info] job_id=${SLURM_JOB_ID}"
echo "[info] array_task_id=${TASK_INDEX}"
echo "[info] variant=${VARIANT_NAME}"
echo "[info] hostname=$(hostname)"
echo "[info] start_time=$(date -Iseconds)"
echo "[info] run_config=${RUN_CONFIG}"
echo "[info] output_root=${OUTPUT_ROOT}"

module purge
module load R/4.4.0-openblas-rocky8 gdal/3.8.2-rocky8 geos/3.11.2-rocky8 udunits/2.2.27.6_gcc11.3.0-rocky8

Rscript "${PIPELINE_ROOT}/Scripts/run_all_eval.R" "${RUN_CONFIG}"

echo "[info] end_time=$(date -Iseconds)"
