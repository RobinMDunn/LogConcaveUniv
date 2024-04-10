# LogConcaveUniv Repository Overview

This repository contains code to replicate the results of [Universal Inference Meets Random Projections: A Scalable Test for Log-concavity](https://arxiv.org/abs/2111.09254) by [Robin Dunn](https://robinmdunn.github.io/), [Aditya Gangrade](https://adityagangrade.wordpress.com/), [Larry Wasserman](https://www.stat.cmu.edu/~larry/), and [Aaditya Ramdas](http://www.stat.cmu.edu/~aramdas/). The `LogConcaveUniv` library contains functions to reproduce the simulations.

## Folder contents

- [R](R): Function definitions
- [man](man): Package documentation files
- [sim_batch_scripts](sim_batch_scripts): Contains SLURM batch scripts to run the simulations. Scripts are labeled by the figure associated with the simulations. These scripts run the code in [sim_code](sim_code), using the parameters in [sim_params](sim_params). 
- [sim_code](sim_code): Code for the paper's simulations. Scripts are labeled by the figure for which they simulate data. Each R script saves the simulation output to [sim_data](sim_data).
- [sim_data](sim_data): Output of simulations from [sim_code](sim_code).
- [sim_params](sim_params): Parameters for simulations. Each row contains a single choice of parameters. The scripts in [sim_code](sim_code) read in these files, and the scripts in [sim_batch_scripts](sim_batch_scripts) loop through all choices of parameters.
- [sim_plot_code](sim_plot_code): Code to reproduce the paper's plots and tables. Reads in data from [sim_data](sim_data) and outputs plots to [sim_plots](sim_plots).
- [sim_plots](sim_plots): Plots from the paper. The plots are the output of the scripts in [sim_plot_code](sim_plot_code).

## Installing the LogConcaveUniv package

We can use `devtools` to install the `LogConcaveUniv` package.

```
# Install and load devtools
install.packages("devtools")
library(devtools)

# Install and load LogConcaveUniv
devtools::install_github("RobinMDunn/LogConcaveUniv")
library(LogConcaveUniv)
```

## Run the simulations for a given figure

### If you have access to a supercomputer with a Slurm workload manager
In the [sim_batch_scripts](sim_batch_scripts) folder, scripts are labeled by their associated figure. Run all batch scripts corresponding to the figure of interest. The allocated run time is estimated from the choice of parameters for which the code has the longest run time. Many scripts will run faster than this time. The files in [sim_code](sim_code) each contain progress bars to estimate the remaining run time. You may wish to start running these files outside of a batch submission to understand the run time on your computing system. 

The simulation output will be stored in the [sim_data](sim_data) folder, with one dataset per choice of parameters. To combine these datasets into a single dataset (as they currently appear in [sim_data](sim_data)), run the code in [sim_code/combine_datasets.R](sim_code/combine_datasets.R).

**Example**: [sim_batch_scripts/fig01_fully_NP_randproj.sh](sim_batch_scripts/fig01_fully_NP_randproj.sh)

This script reproduces the universal test simulations for Figure 1. To do this, it runs the R script at [sim_code/fig01_fully_NP_randproj.R](sim_code/fig01_fully_NP_randproj.R). It reads in the parameters from [sim_params/fig01_fully_NP_randproj_params.csv](sim_params/fig01_fully_NP_randproj_params.csv). There are 30 sets of parameters in total. The results will be stored in the [sim_data](sim_data) folder, with names such as fig01_fully_NP_randproj_1.csv, ..., fig01_fully_NP_randproj_30.csv. To combine these files into a single .csv file, run the code at [sim_code/combine_datasets.R](sim_code/combine_datasets.R).

### Without a supercomputer

To run the code without using a job submission system, click on any .sh file. The Rscript lines can be run on a terminal, replacing $SLURM_ARRAY_TASK_ID with the individual indices in the batch array. 

## Reproduce a figure without rerunning the simulations

The R scripts in [sim_plot_code](sim_plot_code) read in the necessary simulated data from the [sim_data](sim_data) folder and output the figures to the [sim_plots](sim_plots) folder. The name of each R script corresponds to the figure that the script produces.
