#!/bin/bash
#SBATCH --job-name=mainv2_xenium_gene_rate
#SBATCH -p msilarge
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --nice=50
#SBATCH --chdir=/users/9/feng0421/Code/R/Evaluation/MainV2

set -euo pipefail

PIPELINE_ROOT="/users/9/feng0421/Code/R/Evaluation/MainV2"
RUN_CONFIG="${PIPELINE_ROOT}/Config/Runs/xenium_lung_ffpe_shared_gene_rate.R"
OUTPUT_ROOT="/projects/standard/zhangl4/shared/Data_feng0421/Evaluation_MainV2/xenium_lung_ffpe_gene_rate"
LOG_DIR="${OUTPUT_ROOT}/logs"

mkdir -p "${OUTPUT_ROOT}" "${LOG_DIR}"

exec > >(tee -a "${LOG_DIR}/mainv2_xenium_gene_rate_${SLURM_JOB_ID}.log")
exec 2>&1

echo "[info] job_id=${SLURM_JOB_ID}"
echo "[info] hostname=$(hostname)"
echo "[info] start_time=$(date -Iseconds)"
echo "[info] run_config=${RUN_CONFIG}"
echo "[info] output_root=${OUTPUT_ROOT}"

module purge
module load R/4.4.0-openblas-rocky8 gdal/3.8.2-rocky8 geos/3.11.2-rocky8 udunits/2.2.27.6_gcc11.3.0-rocky8

Rscript "${PIPELINE_ROOT}/Scripts/run_all_eval.R" "${RUN_CONFIG}"

echo "[info] end_time=$(date -Iseconds)"
