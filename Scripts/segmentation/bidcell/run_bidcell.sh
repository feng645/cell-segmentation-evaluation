#!/bin/bash
#SBATCH --job-name=bidcell_job          # job name
#SBATCH --cpus-per-task=16                # CPU cores per task
#SBATCH -p msigpu                       # GPU partition
#SBATCH --gres=gpu:1                # one A40 GPU
#SBATCH --time=24:00:00                 # wall time limit
#SBATCH --mem=150GB                      # memory
#SBATCH --chdir=/users/9/feng0421/Code/BIDCell  # working directory

set -euo pipefail

module load miniforge
eval "$(conda shell.bash hook)"    # load the shell functions needed by conda activate
conda activate bidcell
# conda install -c anaconda pyyaml
# conda install -c anaconda pydantic
# conda install numpy<2.0.0
# python -m pip install bidcell
# python -m pip install "numpy<2.0.0"
# python -m pip install cellpose==2.2.3

LOG_DIR="/users/9/feng0421/Code/BIDCell/logs"
mkdir -p "$LOG_DIR"

RUN_TS="$(date +%Y%m%d_%H%M%S)"
TIME_LOG="$LOG_DIR/bidcell_time_${RUN_TS}.log"

echo "Job started: $(date)" | tee "$TIME_LOG"
echo "Host: $(hostname)" | tee -a "$TIME_LOG"
echo "Working directory: $(pwd)" | tee -a "$TIME_LOG"
echo "SLURM job ID: ${SLURM_JOB_ID:-N/A}" | tee -a "$TIME_LOG"
echo "Timing command: /usr/bin/time -v python run_bidcell.py" | tee -a "$TIME_LOG"
echo "" | tee -a "$TIME_LOG"

/usr/bin/time -v -o "$TIME_LOG" -a python run_bidcell.py  # run the BIDCell pipeline

echo "" | tee -a "$TIME_LOG"
echo "Job finished: $(date)" | tee -a "$TIME_LOG"
echo "Time log saved to: $TIME_LOG" | tee -a "$TIME_LOG"
