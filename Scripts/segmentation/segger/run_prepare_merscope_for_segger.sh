#!/bin/bash
#SBATCH --job-name=segger_prepare_merscope
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH -p msismall
#SBATCH --mem=32GB
#SBATCH --nice=99
#SBATCH --chdir=/users/9/feng0421/Code/Segger

module purge
module load miniforge gcc clang
eval "$(conda shell.bash hook)"
conda activate segger

echo "which python:"
which python
python --version

mkdir -p /scratch.global/feng0421/output/MerFish/segger_ready

# python prepare_merscope_for_segger.py \
#   --input-dir /users/9/feng0421/BIDCell/data/MerFish \
#   --output-dir /scratch.global/feng0421/output/MerFish/segger_ready \
#   --overwrite \
#   --min-boundary-span 1.0


python prepare_merscope_for_segger.py \
  --input-dir /users/9/feng0421/BIDCell/data/MerFish \
  --transcripts-csv /users/9/feng0421/BIDCell/data/MerFish/roi_transcripts.csv \
  --output-dir /scratch.global/feng0421/output/MerFish/segger_ready \
  --overwrite \
  --min-boundary-span 1.0
