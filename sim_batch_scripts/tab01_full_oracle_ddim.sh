#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-12:1

Rscript sim_code/tab01_full_oracle_ddim.R sim_params/tab01_full_oracle_ddim_params.csv $SLURM_ARRAY_TASK_ID

