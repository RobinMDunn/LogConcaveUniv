#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-12:1

Rscript sim_code/tab01_fully_NP_axis.R sim_params/tab01_fully_NP_axis_params.csv $SLURM_ARRAY_TASK_ID

