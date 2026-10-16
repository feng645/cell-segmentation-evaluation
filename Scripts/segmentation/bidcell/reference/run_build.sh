#!/bin/bash
#SBATCH --job-name=Build_ref  # job name
#SBATCH --cpus-per-task=16
#SBATCH -p msismall
#SBATCH --time=2:00:00
#SBATCH --mem=100GB
#SBATCH --chdir=/users/9/feng0421/Code/BIDCell/sc_ref  # working directory

echo "Started: $(date)"

# Activate the environment
module load miniforge
eval "$(conda shell.bash hook)"
conda activate bidcell

# Run the script
python build_ref.py --mode create \
  --input_file /scratch.global/feng0421/data/ref/HCLA.h5ad \
  --output_file lung_data.csv \
  --ct_map '{
  "AT2": 0,
  "AT1": 1,
  "Smooth muscle": 2,
  "CD8 T cells": 3,
  "CD4 T cells": 4,
  "T cells proliferating": 5,
  "NK cells": 6,
  "B cells": 7,
  "Plasma cells": 8,
  "Classical monocytes": 9,
  "Alveolar macrophages": 10,
  "Monocyte-derived Mph": 11,
  "Fibroblasts": 12,
  "Endothelial cells": 13
}'

echo "Finished: $(date)"