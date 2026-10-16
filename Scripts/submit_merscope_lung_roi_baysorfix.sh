#!/bin/bash
#SBATCH --job-name=eval_baysorfix_merscope
#SBATCH -p msismall
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --nice=50
#SBATCH --chdir=/users/9/feng0421/Code/R/Evaluation/MainV2
#SBATCH --output=/projects/standard/zhangl4/shared/Data_feng0421/Baysor_rerun_unassigned_fix/logs/%x_%j.out

# MainV2 evaluation with the rerun Baysor polygons; other four methods unchanged.
# Submit with --dependency=afterok:<baysor job id>.
set -euo pipefail
echo "Job ${SLURM_JOB_ID} on $(hostname), started $(date)"
module purge
module load R/4.4.0-openblas-rocky8 gdal/3.8.2-rocky8 geos/3.11.2-rocky8 udunits/2.2.27.6_gcc11.3.0-rocky8
Rscript Scripts/run_all_eval.R Config/Runs/merscope_lung_roi_gene_rate_baysorfix.R
echo "Job finished $(date)"
