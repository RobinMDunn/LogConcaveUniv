#!/bin/bash
#SBATCH -t 05:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-96:1

Rscript sim_code/fig05_fully_NP_randproj.R sim_params/fig05_fully_NP_randproj_params.csv $SLURM_ARRAY_TASK_ID

