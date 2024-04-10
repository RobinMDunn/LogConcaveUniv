#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1

Rscript sim_code/fig10_densities.R

