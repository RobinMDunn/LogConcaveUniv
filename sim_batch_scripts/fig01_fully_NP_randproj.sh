#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-55:1

Rscript sim_code/fig01_fully_NP_randproj.R sim_params/fig01_fully_NP_randproj_params.csv $SLURM_ARRAY_TASK_ID

