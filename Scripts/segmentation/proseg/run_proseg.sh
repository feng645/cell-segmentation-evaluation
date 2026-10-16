#!/bin/bash
#SBATCH --job-name=proseg
#SBATCH --cpus-per-task=18
#SBATCH --time=8:00:00
#SBATCH -p msilarge
#SBATCH --mem=150GB
#SBATCH --nice=99
#SBATCH --chdir=/users/9/feng0421/Code/Proseg
set -euo pipefail

# Load modules
module purge
module load miniforge gcc clang
eval "$(conda shell.bash hook)"
conda activate proseg

# Environment variables
export RUST_BACKTRACE=1
export OMP_NUM_THREADS=16
export PROSEG_THREADS=16    

# Paths
INPUT_FILE="/projects/standard/zhangl4/shared/Data_feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts_cleaned.csv"
# INPUT_FILE="/users/9/feng0421/BIDCell/data/CosMx/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_tx_file.csv"
OUTPUT_DIR="/scratch.global/feng0421/output/Xenium_Proseg"
LOG_FILE="proseg_$(date +%Y%m%d_%H%M%S).log"

# Create the output directory
mkdir -p "$OUTPUT_DIR"
LOG_PATH="$OUTPUT_DIR/$LOG_FILE"

format_duration() {
    local total_seconds="$1"
    printf '%02d:%02d:%02d\n' \
        $((total_seconds / 3600)) \
        $(((total_seconds % 3600) / 60)) \
        $((total_seconds % 60))
}

# Job information
echo "======================================"
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "Output directory: $OUTPUT_DIR"
echo "======================================"

# Check the input file
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

echo "Input file: $INPUT_FILE"
echo "File size: $(du -h "$INPUT_FILE" | cut -f1)"
echo "File lines: $(wc -l < "$INPUT_FILE")"

# Check the conda environment and the Proseg installation
echo "Checking environment..."
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Proseg version: $(proseg --version 2>/dev/null || echo "Version info not available")"

# Check the input file format
echo "Checking input file format..."
echo "Header: $(head -1 "$INPUT_FILE")"

# Run Proseg
echo "======================================"
echo "Starting Proseg analysis at: $(date)"
echo "======================================"

# Time the run and send all output to the log file
# time proseg \
#     --cosmx "$INPUT_FILE" \
#     --x-column x_global_px \
#     --y-column y_global_px \
#     --gene-column target \
#     --cell-id-column cell_ID \
#     --ncomponents 7 \
#     --nthreads 16 \
#     --output-expected-counts expected-counts.mtx.gz \
#     --output-cell-metadata cell-metadata.csv.gz \
#     --output-cell-polygons cell-polygons.geojson.gz \
#     --output-path "$OUTPUT_DIR" \
#     2>&1 | tee "$OUTPUT_DIR/$LOG_FILE"


# time proseg \
#     --merscope "$INPUT_FILE" \
#     --ncomponents 7 \
#     --nthreads 16 \
#     --voxel-layers 7 \
#     --output-expected-counts expected-counts.mtx.gz \
#     --output-cell-metadata cell-metadata.csv.gz \
#     --output-cell-polygons cell-polygons.geojson.gz \
#     --output-path "$OUTPUT_DIR" \
#     2>&1 | tee "$OUTPUT_DIR/$LOG_FILE"


{
    echo "Segmentation step started: $(date)"
    echo "Host: $(hostname)"
    echo "Working directory: $(pwd)"
    echo "SLURM job ID: ${SLURM_JOB_ID:-N/A}"
    echo "Timing command: /usr/bin/time -v proseg ..."
    echo
} > "$LOG_PATH"

start_epoch=$(date +%s)
proseg_status=0
if ! /usr/bin/time -v proseg \
    --xenium "$INPUT_FILE" \
    --ncomponents 7 \
    --nthreads 16 \
    --output-expected-counts expected-counts.mtx.gz \
    --output-cell-metadata cell-metadata.csv.gz \
    --output-cell-polygons cell-polygons.geojson.gz \
    --output-path "$OUTPUT_DIR" \
    2>&1 | tee -a "$LOG_PATH"; then
    proseg_status=${PIPESTATUS[0]}
fi
end_epoch=$(date +%s)
runtime_seconds=$((end_epoch - start_epoch))
runtime_hms=$(format_duration "$runtime_seconds")
peak_mem_kb=$(grep "Maximum resident set size" "$LOG_PATH" | awk -F: '{gsub(/^[ \t]+/, "", $2); print $2}' | tail -n1)

{
    echo
    echo "======================================"
    echo "Proseg analysis completed at: $(date)"
    echo "Output directory: $OUTPUT_DIR"
    echo "Runtime (seconds): $runtime_seconds"
    echo "Runtime (HH:MM:SS): $runtime_hms"
    echo "Peak memory (kbytes): ${peak_mem_kb:-N/A}"
    echo "Exit status: $proseg_status"
} | tee -a "$LOG_PATH"

exit "$proseg_status"
