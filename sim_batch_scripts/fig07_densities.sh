#!/bin/bash
#SBATCH -t 00:01:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1

Rscript sim_code/fig07_densities.R

