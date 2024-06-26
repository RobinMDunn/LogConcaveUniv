# Test H_0: log-concave versus H_1: not log-concave using universal LRT.
# Split data into D_0 and D_1.
# Use true density as numerator. Get log-concave MLE on D_0.
# Evaluate at all values in D_0.
# Density is 0.5*N(0, I_2) + 0.5*N(0, sigma^2 I_2), where sigma = sqrt(3).

# Read in library
library(LogConcaveUniv)

# Read in arguments for file with all parameters,
# line number for parameters for current simulation, and
# number of parallel cores to use in simulation.

parameter_file <- "sim_params/fig06_full_oracle_ddim_params.csv"
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
                                  sigma = sigma, B = B,
                                  sim = start_sim:(start_sim + n_sim - 1),
                                  alpha = alpha, p_0 = p,
                                  compute_ts = compute_ts)

# Code to run one simulation to check whether to reject H_0
one_sim_full_oracle_ddim <- function(n_obs, d, sigma, B, sim, 
                                     alpha, p_0, compute_ts, mu) {

  # Generate sample from two-component normal mixture model
  true_sample <- matrix(NA, nrow = n_obs, ncol = d)

  for(i in 1:n_obs) {
    mixture_comp <- rbinom(n = 1, size = 1, prob = p_0)
    if(mixture_comp == 0) {
      true_sample[i, ] <- rnorm(n = d, mean = 0, sd = 1)
    } else if(mixture_comp == 1) {
      true_sample[i, ] <- rnorm(n = d, mean = 0, sd = sigma)
    }
  }

  # Run full oracle d-dimensional test to determine whether to reject H_0
  test_out <- LogConcaveUniv::full_oracle_ddim(data = true_sample, B = B,
                                               alpha = alpha, mu = mu,
                                               sigma = sigma, p = p_0,
                                               compute_ts = compute_ts)
  
  return(test_out)

}

# Run simulations to check whether to reject H_0, iterating over rows of results
test_out <- clustermq::Q_rows(df = results, fun = one_sim_full_oracle_ddim, 
                              n_jobs = n_cores, const = list(mu = mu))

# Append outputs to results df
results$avg_ts <- sapply(test_out, FUN = function(x) x$test_stat)
results$reject <- sapply(test_out, FUN = function(x) x$reject_null)

# Create aggregated results (number of rejections at given parameters)
results_agg <- results %>% 
  dplyr::group_by(n_obs, d, sigma, B, alpha, p_0) %>% 
  dplyr::summarize(n_reject = sum(reject), 
                   n_sim = dplyr::n())

# Save simulation results
data.table::fwrite(results_agg, 
                   file = paste0("sim_data/fig06_full_oracle_ddim_",
                                 line_number, ".csv"))
