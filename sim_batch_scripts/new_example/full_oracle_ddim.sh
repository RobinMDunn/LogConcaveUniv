#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-10:1

Rscript sim_code/new_example/sim_full_oracle_ddim.R sim_params/new_example/full_oracle_ddim_params.csv $SLURM_ARRAY_TASK_ID 1

