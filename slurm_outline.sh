#!/bin/bash

#SBATCH --job-name=sims_08022026                                                # Name of the job
#SBATCH --mail-user=agfoes@unc.edu                                              # Email address
#SBATCH --mail-type=ALL                                                         # Alerts sent when job begins, ends, or aborts
#SBATCH --time=48:00:00                                                         # Maximum run time
#SBATCH --nodes=1                                                               # Number of nodes
#SBATCH --ntasks-per-node=1                                                     # Number of tasks per node
#SBATCH --cpus-per-task=1                                                       # Number of CPU cores per task
#SBATCH --mem=24G                                                               # Memory per node
#SBATCH --array=1-500%100                                                           # Array range for jobs
#SBATCH --output=/work/users/a/g/agfoes/P1/results/expanded_rep_sims_08022026/logs/%A_%a.out   # Standard output log
#SBATCH --error=/work/users/a/g/agfoes/P1/results/expanded_rep_sims_08022026/logs/%A_%a.err    # Standard error log

## add R module
module add r/4.5.0

R CMD BATCH --no-restore \
  /work/users/a/g/agfoes/P1/run_simulation.R \
  /work/users/a/g/agfoes/P1/results/expanded_rep_sims_08022026/Rout/run_simulation_$SLURM_ARRAY_TASK_ID.Rout
