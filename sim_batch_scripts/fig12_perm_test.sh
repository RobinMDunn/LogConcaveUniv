#!/bin/bash
#SBATCH -t 17:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-150:1

Rscript sim_code/fig12_perm_test.R sim_params/fig12_perm_test_params.csv $SLURM_ARRAY_TASK_ID

