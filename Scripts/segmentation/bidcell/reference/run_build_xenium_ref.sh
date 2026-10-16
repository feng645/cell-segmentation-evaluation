#!/bin/bash
#SBATCH --job-name=xenium_ref
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH -p msismall
#SBATCH --mem=64GB
#SBATCH --nice=99
#SBATCH --chdir=/users/9/feng0421/Code/BIDCell/sc_ref

set -euo pipefail

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"

module load miniforge
eval "$(conda shell.bash hook)"
conda activate bidcell

TRANSCRIPTS="/scratch.global/feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts.csv.gz"
REF_DIR="/users/9/feng0421/BIDCell/data/sc_references"
SCRIPT="/users/9/feng0421/Code/BIDCell/sc_ref/build_xenium_reference.py"

python "$SCRIPT" \
  --transcripts "$TRANSCRIPTS" \
  --hcla "$REF_DIR/HCLA.h5ad" \
  --nsclc \
    "$REF_DIR/NSCLC_GSE148071_expression_Celltype_majorlineage.txt" \
  --output-dir "$REF_DIR"

# Optional: clean Xenium transcripts using the generated filter list.
python /users/9/feng0421/Code/BIDCell/sc_ref/clean_trans.py \
  "$TRANSCRIPTS" \
  /scratch.global/feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts_filtered.csv \
  -f "$REF_DIR/xenium_target_to_filter.txt" \
  --target-column feature_name

echo "Job finished: $(date)"
