#!/bin/bash
#SBATCH -t 11:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-84:1

Rscript sim_code/fig05_partial_oracle_randproj.R sim_params/fig05_partial_oracle_randproj_params.csv $SLURM_ARRAY_TASK_ID

