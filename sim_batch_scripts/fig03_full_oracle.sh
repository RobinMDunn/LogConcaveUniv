#!/bin/bash
#SBATCH -t 17:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-84:1

Rscript sim_code/fig03_full_oracle.R sim_params/fig03_full_oracle_params.csv $SLURM_ARRAY_TASK_ID

