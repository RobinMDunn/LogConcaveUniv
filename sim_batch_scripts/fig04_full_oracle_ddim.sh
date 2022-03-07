#!/bin/bash
#SBATCH -t 17:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-106:1

Rscript sim_code/fig04_full_oracle_ddim.R sim_params/fig04_full_oracle_ddim_params.csv $SLURM_ARRAY_TASK_ID

