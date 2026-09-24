#!/bin/bash

#SBATCH --job-name=multimodal_testing_P1
#SBATCH --mail-user=agfoes@unc.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-300
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/work/users/a/g/agfoes/P1/results/multimodal_tests/logs/gelc_%A_%a.out
#SBATCH --error=/work/users/a/g/agfoes/P1/results/multimodal_tests/logs/gelc_%A_%a.err

module purge
module load r/4.5.0

Rscript \
  /work/users/a/g/agfoes/P1/run_R/run_multimodal_sim_P1.R