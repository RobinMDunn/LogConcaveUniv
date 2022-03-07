#!/bin/bash
#SBATCH -t 05:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-55:1

Rscript sim_code/fig01_perm_test.R sim_params/fig01_perm_test_params.csv $SLURM_ARRAY_TASK_ID

