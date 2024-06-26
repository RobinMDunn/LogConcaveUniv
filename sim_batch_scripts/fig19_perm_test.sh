#!/bin/bash
#SBATCH -t 04:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-60:1

Rscript sim_code/fig19_perm_test.R sim_params/fig19_perm_test_params.csv $SLURM_ARRAY_TASK_ID

