#!/bin/bash
#SBATCH --job-name=clean_xenium_tx
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32GB
#SBATCH -p msismall
#SBATCH --nice=50
#SBATCH --chdir=/users/9/feng0421/Code/Baysor

set -euo pipefail

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "SLURM job ID: ${SLURM_JOB_ID:-N/A}"

module load miniforge
eval "$(conda shell.bash hook)"
conda activate bidcell

INPUT_TX="/scratch.global/feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts.csv.gz"
OUTPUT_TX="/projects/standard/zhangl4/shared/Data_feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts_cleaned.csv"
FILTER_TXT="/users/9/feng0421/BIDCell/data/sc_references/xenium_target_to_filter.txt"
SCRIPT="/users/9/feng0421/Code/BIDCell/sc_ref/clean_trans.py"

echo "Input transcripts: $INPUT_TX"
echo "Output transcripts: $OUTPUT_TX"
echo "Generated filter list: $FILTER_TXT"
echo "Filtering prefixes: NegControlProbe_ NegControlCodeword_ UnassignedCodeword_ antisense_ BLANK_ Blank- NegPrb"

python "$SCRIPT" \
  "$INPUT_TX" \
  "$OUTPUT_TX" \
  --target-column feature_name \
  --prefixes NegControlProbe_ NegControlCodeword_ UnassignedCodeword_ antisense_ BLANK_ Blank- NegPrb \
  --generated-targets-file "$FILTER_TXT"

echo "Job finished: $(date)"
