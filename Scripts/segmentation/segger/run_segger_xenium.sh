#!/bin/bash
#SBATCH --job-name=segger_xenium
#SBATCH --cpus-per-task=18
#SBATCH --time=12:00:00
#SBATCH -p msigpu
#SBATCH --mem=180GB
#SBATCH --nice=99
#SBATCH --gres=gpu:1
#SBATCH --chdir=/users/9/feng0421/Code/Segger

module purge
module load miniforge gcc clang
eval "$(conda shell.bash hook)"
conda activate segger

# set -euo pipefail

# Ensure CuPy is installed for the CUDA version on this GPU node.
if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "nvidia-smi not found. Please run this job on a GPU node with NVIDIA drivers." >&2
  exit 1
fi

CUDA_MAJOR=$(awk -F'CUDA Version: ' 'NF>1 {split($2,a," "); split(a[1],b,"."); print b[1]; exit}' <(nvidia-smi))
case "${CUDA_MAJOR}" in
  13) CUPY_PKG="cupy-cuda12x" ;;
  12) CUPY_PKG="cupy-cuda12x" ;;
  11) CUPY_PKG="cupy-cuda11x" ;;
  *)
    echo "Unsupported CUDA major version from nvidia-smi: ${CUDA_MAJOR}" >&2
    exit 1
    ;;
esac

if pip show "${CUPY_PKG}" >/dev/null 2>&1; then
  echo "CuPy (${CUPY_PKG}) already installed in environment segger."
else
  echo "Installing ${CUPY_PKG} into environment segger..."
  pip install "${CUPY_PKG}<14" "numpy>=1.26,<2"
fi

echo "Installing CUDA runtime/NVRTC Python wheels..."
pip install -U "nvidia-cuda-nvrtc-cu12>=12,<13" "nvidia-cuda-runtime-cu12>=12,<13"

NVIDIA_LIBS=$(
python - <<'PY'
import glob
import site

paths = []
for p in site.getsitepackages():
    paths.extend(glob.glob(p + "/nvidia/cuda_nvrtc/lib"))
    paths.extend(glob.glob(p + "/nvidia/cuda_runtime/lib"))
print(":".join(paths))
PY
)
if [[ -n "${NVIDIA_LIBS}" ]]; then
  export LD_LIBRARY_PATH="${NVIDIA_LIBS}:${LD_LIBRARY_PATH}"
fi

python - <<'PY'
import ctypes
ctypes.CDLL("libnvrtc.so.12")
print("OK: libnvrtc.so.12 found")
PY

python - <<'PY'
import cupy as cp
print("CuPy:", cp.__version__)
print("Visible GPU count:", cp.cuda.runtime.getDeviceCount())
PY

export TMPDIR=/tmp
export NUMBA_CACHE_DIR=/tmp/numba_cache_${USER}
mkdir -p "$NUMBA_CACHE_DIR"

SEGGER_CLI_ROOT=/users/9/feng0421/segger_dev/src/segger/cli
XENIUM_BASE=/scratch.global/feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs
SEGGER_OUT_ROOT=/scratch.global/feng0421/output/Xenium_segger
DATA_TIDY_DIR=${SEGGER_OUT_ROOT}/data_tidy
MODELS_DIR=${SEGGER_OUT_ROOT}/models
PRED_DIR=${SEGGER_OUT_ROOT}/predictions
LOG_DIR=${SEGGER_OUT_ROOT}/logs
PEAK_RAM_FILE=${LOG_DIR}/peak_ram_summary.txt

mkdir -p "$DATA_TIDY_DIR" "$MODELS_DIR" "$PRED_DIR" "$LOG_DIR"
: > "$PEAK_RAM_FILE"

run_step() {
  local step_name="$1"
  shift
  local log_file="${LOG_DIR}/${step_name}.log"

  echo "[$(date '+%F %T')] START ${step_name}" | tee -a "${LOG_DIR}/run_summary.log"
  set +e
  /usr/bin/time -v "$@" &> "$log_file"
  local rc=$?
  set -e

  local peak_kb
  peak_kb=$(awk -F: '/Maximum resident set size/ {gsub(/^[ \t]+/, "", $2); print $2}' "$log_file" | tail -n 1)
  if [[ -z "${peak_kb}" ]]; then
    peak_kb="NA"
  fi
  echo "${step_name} peak_rss_kb=${peak_kb}" | tee -a "$PEAK_RAM_FILE"
  echo "[$(date '+%F %T')] END ${step_name} rc=${rc}" | tee -a "${LOG_DIR}/run_summary.log"

  if [[ $rc -ne 0 ]]; then
    echo "Step failed: ${step_name}. See ${log_file}" >&2
    exit $rc
  fi
}

# run_step step1_create_dataset \
#   python ${SEGGER_CLI_ROOT}/create_dataset_fast.py \
#   --base_dir ${XENIUM_BASE} \
#   --data_dir ${DATA_TIDY_DIR} \
#   --sample_type xenium \
#   --scrnaseq_file /scratch.global/feng0421/output/MerFish/segger_ready/HCLA_symbol_varnames.h5ad \
#   --celltype_column ann_level_3 \
#   --k_bd 4 \
#   --dist_bd 15 \
#   --k_tx 3 \
#   --dist_tx 5 \
#   --tile_size 10000 \
#   --neg_sampling_ratio 5 \
#   --frac 1 \
#   --val_prob 0.1 \
#   --test_prob 0.2 \
#   --n_workers 16

# run_step step2_train_model \
#   python ${SEGGER_CLI_ROOT}/train_model.py \
#   --dataset_dir ${DATA_TIDY_DIR} \
#   --models_dir ${MODELS_DIR} \
#   --sample_tag lung_cancer_v1 \
#   --init_emb 8 \
#   --hidden_channels 32 \
#   --num_mid_layers 2 \
#   --out_channels 8 \
#   --heads 2 \
#   --batch_size 4 \
#   --num_workers 4 \
#   --accelerator cuda \
#   --devices 1 \
#   --strategy auto \
#   --precision 16-mixed \
#   --max_epochs 200 \
#   --learning_rate 1e-3 \
#   --save_best_model True

run_step step3_predict \
  python ${SEGGER_CLI_ROOT}/predict_fast.py \
  --segger_data_dir ${DATA_TIDY_DIR} \
  --models_dir ${MODELS_DIR} \
  --benchmarks_dir ${PRED_DIR} \
  --transcripts_file ${XENIUM_BASE}/transcripts.parquet \
  --save_tag lung_cancer_v1_pred \
  --model_version 0 \
  --batch_size 8 \
  --num_workers 2 \
  --use_cc True \
  --min_transcripts 5 \
  --cell_id_col segger_cell_id \
  --knn_method kd_tree \
  --file_format anndata \
  --k_bd 4 \
  --dist_bd 6 \
  --k_tx 5 \
  --dist_tx 5

echo "All steps completed. Peak RAM summary: ${PEAK_RAM_FILE}"
