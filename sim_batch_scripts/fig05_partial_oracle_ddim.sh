#!/bin/bash
#SBATCH -t 23:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-65:1

Rscript sim_code/fig05_partial_oracle_ddim.R sim_params/fig05_partial_oracle_ddim_params.csv $SLURM_ARRAY_TASK_ID

