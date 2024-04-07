#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-8:1

Rscript sim_code/fig13_test_stats.R sim_params/fig13_test_stats_params.csv $SLURM_ARRAY_TASK_ID

