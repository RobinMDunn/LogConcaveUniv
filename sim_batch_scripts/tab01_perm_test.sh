#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-12:1

Rscript sim_code/tab01_perm_test.R sim_params/tab01_perm_test_params.csv $SLURM_ARRAY_TASK_ID

