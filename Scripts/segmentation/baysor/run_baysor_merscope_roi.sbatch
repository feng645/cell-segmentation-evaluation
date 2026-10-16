#!/bin/bash
#SBATCH --job-name=baysor_fix_merscope
#SBATCH -p msismall
#SBATCH --constraint=genoa  # libbaysor.so sysimage was built on genoa (original runs: n54, n40)
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --nice=5
#SBATCH --chdir=/users/9/feng0421/Code/Baysor
#SBATCH --output=/projects/standard/zhangl4/shared/Data_feng0421/Baysor_rerun_unassigned_fix/logs/%x_%j.out

# Rerun of Baysor on the MERSCOPE ROI with the unassigned prior label fixed ("-1").
# All other parameters are identical to the original run
# (Code/R/Evaluation/Baysor_out/MERFISH_Baysor_params.dump.toml).
set -euo pipefail

ROOT=/projects/standard/zhangl4/shared/Data_feng0421/Baysor_rerun_unassigned_fix
INPUT=/users/9/feng0421/BIDCell/data/MerFish/roi_transcripts.csv
OUT_DIR="$ROOT/MERSCOPE_ROI"
CONFIG="$ROOT/configs/merscope_roi.toml"
TIME_LOG="$ROOT/logs/baysor_merscope_time_${SLURM_JOB_ID}.log"

export PATH="$HOME/.juliaup/bin:$PATH"
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Job ${SLURM_JOB_ID} on $(hostname), started $(date)"
mkdir -p "$OUT_DIR"

/usr/bin/time -v -o "$TIME_LOG" ./baysor run \
    "$INPUT" \
    :cell_id \
    --config "$CONFIG" \
    --output "$OUT_DIR/segmentation.csv" \
    --x-column global_x \
    --y-column global_y \
    --z-column global_z \
    --gene-column gene \
    --min-molecules-per-cell 10 \
    --n-clusters 7 \
    --count-matrix-format "tsv" \
    --prior-segmentation-confidence 0.5 \
    --plot

# Fail loudly if Baysor still did not find any unassigned molecules.
if grep -q "No unassigned molecules found" "$OUT_DIR/segmentation_log.log"; then
    echo "ERROR: unassigned prior label still not recognized" >&2
    exit 1
fi
echo "Job finished $(date)"
