# Test H_0: log-concave versus H_1: not log-concave using universal LRT.
# Choose a random point along surface of d-dimensional unit sphere.
# Project observations onto this vector, and rotate to x-axis.
# Check whether this one-dimensional projection is log-concave.
# Split data into D_0 and D_1.
# Get kernel density estimate on D_1. Get log-concave MLE on D_0.
# Evaluate at all values in D_0.
# Repeat at multiple projections. Reject if mean test stat >= 1/alpha.
# Density is 0.5*N(0, I_2) + 0.5*N(0, sigma^2 I_2), where sigma = sqrt(3).

# Read in library
library(LogConcaveUniv)

# Read in arguments for file with all parameters,
# line number for parameters for current simulation, and
# number of parallel cores to use in simulation.

parameter_file <- "sim_params/new_example/fully_NP_randproj_params.csv"
line_number <- 1
n_cores <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  parameter_file <- args[1]
  line_number <- as.numeric(args[2])
  n_cores <- as.numeric(args[3])
}

parameter_df <- data.table::fread(parameter_file)

# Assign arguments based on input.
# Arguments are d (dimension),
# sigma (generating from (1-p) N(0, I_d) + p N(0, sigma^2 I_d)),
# n_obs (number of obs),
# start_sim (index of starting sim), n_sim (number of sims),
# n_proj (number of random projections), 
# B (number of subsamples),
# compute_ts (indicator for whether to compute test stat.
#             If 1, computes test stat.
#             If 0, does not compute test stat, allowing for early rejection.)

parameters <- parameter_df %>% dplyr::slice(line_number)
d <- parameters$d
sigma <- parameters$sigma
n_obs <- parameters$n_obs
start_sim <- parameters$start_sim
n_sim <- parameters$n_sim
n_proj <- parameters$n_proj
B <- parameters$B
compute_ts <- parameters$compute_ts

# Create mu vector. (mu is the 0 vector for both components.)
mu <- rep(0, d)

# Proportions for the two Gaussian components
p <- 0.5

# Alpha level
alpha <- 0.1

# Create data frame to store results
results <- data.table::data.table(n_obs = n_obs, d = d, 
                                  sigma = sigma,
                                  n_proj = n_proj, B = B,
                                  sim = start_sim:(start_sim + n_sim - 1),
                                  alpha = alpha, p_0 = p,
                                  compute_ts = compute_ts)

# Code to run one simulation to check whether to reject H_0
one_sim_fully_NP_randproj <- function(n_obs, d, sigma, n_proj, B, sim, 
                                      alpha, p_0, compute_ts) {
  
  # Generate sample from two-component normal location model
  true_sample <- matrix(NA, nrow = n_obs, ncol = d)

  for(i in 1:n_obs) {
    mixture_comp <- rbinom(n = 1, size = 1, prob = p_0)
    if(mixture_comp == 0) {
      true_sample[i, ] <- rnorm(n = d, mean = 0, sd = 1)
    } else if(mixture_comp == 1) {
      true_sample[i, ] <- rnorm(n = d, mean = 0, sd = sigma)
    }
  }

  # Run fully nonparametric random projection test to determine
  # whether to reject H_0
  test_out <- LogConcaveUniv::fully_NP_randproj(data = true_sample, B = B,
                                                n_proj = n_proj, alpha = alpha,
                                                compute_ts = compute_ts)

  return(test_out)

}

# Run simulations to check whether to reject H_0, iterating over rows of results
test_out <- clustermq::Q_rows(df = results, fun = one_sim_fully_NP_randproj, 
                              job_size = n_cores)

# Append outputs to results df
results$avg_ts <- sapply(test_out, FUN = function(x) x$test_stat)
results$reject <- sapply(test_out, FUN = function(x) x$reject_null)

# Create aggregated results (number of rejections at given parameters)
results_agg <- results %>% 
  dplyr::group_by(n_obs, d, sigma, n_proj, B, alpha, p_0) %>% 
  dplyr::summarize(n_reject = sum(reject), 
                   n_sim = dplyr::n())

# Save simulation results
data.table::fwrite(results_agg, 
                   file = paste0("sim_data/new_example/fully_NP_randproj_",
                                 line_number, ".csv"))
