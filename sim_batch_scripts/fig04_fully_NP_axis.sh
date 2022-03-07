#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-92:1

Rscript sim_code/fig04_fully_NP_axis.R sim_params/fig04_fully_NP_axis_params.csv $SLURM_ARRAY_TASK_ID

