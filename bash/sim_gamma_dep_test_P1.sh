#!/bin/bash

#SBATCH --job-name=centered_X3
#SBATCH --mail-user=agfoes@unc.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-2000%1000
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/work/users/a/g/agfoes/P1/results/dependence_gelc/logs/p1_%A_%a.out
#SBATCH --error=/work/users/a/g/agfoes/P1/results/dependence_gelc/logs/p1_%A_%a.err

module purge
module load r/4.5.0

Rscript \
  /work/users/a/g/agfoes/P1/sim_gamma_dependence_testing_P1.R