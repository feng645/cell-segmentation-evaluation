#!/bin/bash
#SBATCH --job-name=baysor_fix_xenium
#SBATCH -p msismall
#SBATCH --constraint=genoa  # libbaysor.so sysimage was built on genoa (original runs: n54, n40)
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --time=12:00:00
#SBATCH --nice=5
#SBATCH --chdir=/users/9/feng0421/Code/Baysor
#SBATCH --output=/projects/standard/zhangl4/shared/Data_feng0421/Baysor_rerun_unassigned_fix/logs/%x_%j.out

# Rerun of Baysor on the full Xenium section with the unassigned prior label fixed
# ("UNASSIGNED"). All other parameters are identical to the original run
# (Data_feng0421/Xenium_Baysor/segmentation_params.dump.toml).
# Original run: ~1.6-2.5 h wall time, ~45-48 GB peak RSS.
set -euo pipefail

ROOT=/projects/standard/zhangl4/shared/Data_feng0421/Baysor_rerun_unassigned_fix
INPUT=/projects/standard/zhangl4/shared/Data_feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts_cleaned.csv
OUT_DIR="$ROOT/Xenium"
CONFIG="$ROOT/configs/xenium.toml"
TIME_LOG="$ROOT/logs/baysor_xenium_time_${SLURM_JOB_ID}.log"

export PATH="$HOME/.juliaup/bin:$PATH"
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Job ${SLURM_JOB_ID} on $(hostname), started $(date)"
mkdir -p "$OUT_DIR"

/usr/bin/time -v -o "$TIME_LOG" ./baysor run \
    "$INPUT" \
    :cell_id \
    --config "$CONFIG" \
    --output "$OUT_DIR/segmentation.csv" \
    --x-column x_location \
    --y-column y_location \
    --z-column z_location \
    --gene-column feature_name \
    --min-molecules-per-cell 10 \
    --n-clusters 7 \
    --count-matrix-format "tsv" \
    --prior-segmentation-confidence 0.5 \
    --plot

if grep -q "No unassigned molecules found" "$OUT_DIR/segmentation_log.log"; then
    echo "ERROR: unassigned prior label still not recognized" >&2
    exit 1
fi
echo "Job finished $(date)"
