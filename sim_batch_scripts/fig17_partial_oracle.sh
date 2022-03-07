#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --array=1-60:1

Rscript sim_code/fig17_partial_oracle.R sim_params/fig17_partial_oracle_params.csv $SLURM_ARRAY_TASK_ID

