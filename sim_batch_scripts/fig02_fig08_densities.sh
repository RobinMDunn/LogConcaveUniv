#!/bin/bash
#SBATCH -t 04:00:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1

Rscript sim_code/fig02_fig08_densities.R

