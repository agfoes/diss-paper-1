#!/bin/bash

#SBATCH --job-name=gelc_stress_test
#SBATCH --mail-user=agfoes@unc.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-1000%500
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/work/users/a/g/agfoes/P1/results/gelc_comparison_09022026/logs/gelc_%A_%a.out
#SBATCH --error=/work/users/a/g/agfoes/P1/results/gelc_comparison_09022026/logs/gelc_%A_%a.err

module purge
module load r/4.5.0

Rscript \
  /work/users/a/g/agfoes/P1/R/GELc/GELc_testing.R