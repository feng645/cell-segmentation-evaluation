#!/bin/bash
#SBATCH --job-name=find_roi  # job name
#SBATCH --cpus-per-task=2                  
#SBATCH --time=2:00:00
#SBATCH --mem=64GB
#SBATCH --chdir=/users/9/feng0421/Code/Cellpose  # working directory

echo "Started: $(date)"

# Activate the environment
module load miniforge
eval "$(conda shell.bash hook)"
conda activate cellpose

# Run the script
python roi_crop_merfish.py \
  --transcripts /scratch.global/feng0421/data/MERFISH/detected_transcripts_seg_7.csv \
  --dapi /scratch.global/feng0421/data/MERFISH/HumanLungCancerPatient1_images_mosaic_DAPI_z0.tif \
  --transform /scratch.global/feng0421/data/MERFISH/HumanLungCancerPatient1_images_micron_to_mosaic_pixel_transform.csv \
  --roi-um 550 550 800 800 \
  --out-img /users/9/feng0421/Code/R/Evaluation/Cellpose_out/roi_dapi.tif \
  --out-dir /users/9/feng0421/Code/R/Evaluation/Cellpose_out/ \
  --viz True

echo "Finished: $(date)"
