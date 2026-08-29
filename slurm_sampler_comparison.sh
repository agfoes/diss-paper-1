#!/bin/bash

#SBATCH --job-name=sampler_test_custom_joint                                                # Name of the job
#SBATCH --mail-user=agfoes@unc.edu                                              # Email address
#SBATCH --mail-type=ALL                                                         # Alerts sent when job begins, ends, or aborts
#SBATCH --time=48:00:00                                                         # Maximum run time
#SBATCH --nodes=1                                                               # Number of nodes
#SBATCH --ntasks-per-node=1                                                     # Number of tasks per node
#SBATCH --cpus-per-task=1                                                       # Number of CPU cores per task
#SBATCH --mem=16G                                                               # Memory per node
#SBATCH --array=1-100                                                           # Array range for jobs
#SBATCH --output=/work/users/a/g/agfoes/P1/results/sampler_test_custom_joint/logs/%A_%a.out   # Standard output log
#SBATCH --error=/work/users/a/g/agfoes/P1/results/sampler_test_custom_joint/logs/%A_%a.err    # Standard error log

## add R module
module add r/4.5.0

R CMD BATCH --no-restore /work/users/a/g/agfoes/P1/run_sampler_comparison.R /work/users/a/g/agfoes/P1/results/sampler_test_custom_joint/Rout/run_sampler_comparison_$SLURM_ARRAY_TASK_ID.Rout
