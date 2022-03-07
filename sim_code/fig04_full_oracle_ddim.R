# Test H_0: log-concave versus H_1: not log-concave using universal LRT.
# Split data into D_0 and D_1. 
# Use true density as numerator. Get log-concave MLE on D_0.
# Evaluate at all values in D_0.

suppressMessages(library(LogConcDEAD))
suppressMessages(library(MASS))
suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(tidyverse))
suppressMessages(library(mvtnorm))

# Read in arguments for file with all parameters and 
# line number for parameters for current simulation.

parameter_file <- "sim_params/fig04_full_oracle_ddim_params.csv"
line_number <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  parameter_file <- args[1]
  line_number <- as.numeric(args[2])
}

parameter_df <- fread(parameter_file)

# Assign arguments based on input.
# Arguments are d (dimension), mu_norm, n_obs (number of obs),
# start_sim (index of starting sim), n_sim (number of sims), 
# B (number of subsamples),
# equal_space_mu (indicator for mu structure in second component. 
#                 If 1, each component is ||mu||*d^(-1/2).
#                 If 0, mu = (||mu||, 0, ... 0),
# compute_ts (indicator for whether to compute test stat.
#             If 1, computes test stat.
#             If 0, does not compute test stat, allowing for early rejection.)

parameters <- parameter_df %>% slice(line_number)
d <- parameters$d
mu_norm <- parameters$mu_norm
n_obs <- parameters$n_obs
start_sim <- parameters$start_sim
n_sim <- parameters$n_sim
B <- parameters$B
equal_space_mu <- parameters$equal_space_mu
compute_ts <- parameters$compute_ts

# Create mu vector
if(equal_space_mu == 0) {
  mu <- rep(0, d)
  mu[1] <- mu_norm
} else if(equal_space_mu == 1) {
  mu <- rep(mu_norm * d^(-1/2), d)
}

# Proportions for the two Gaussian components
p <- 0.5

# Alpha level
alpha <- 0.1

# Create data frame to store results
results <- data.table(n_obs = n_obs, d = d, mu_norm = mu_norm, 
                      equal_space_mu = equal_space_mu, B = B,
                      sim = start_sim:(start_sim + n_sim - 1), alpha = alpha,
                      p_0 = p, test_stat = NA_real_, reject = NA_real_)

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results), 
                       clear = T, show_after = 0)

# Run simulations to check whether to reject H_0
for(row in 1:nrow(results)) {
  
  # Increment progress bar
  pb$tick()
  
  # Generate sample from two-component normal location model
  true_sample <- matrix(NA, nrow = n_obs, ncol = d)
  
  for(i in 1:n_obs) {
    mixture_comp <- rbinom(n = 1, size = 1, prob = p)
    if(mixture_comp == 0) {
      true_sample[i, ] <- rnorm(n = d, mean = 0, sd = 1)
    } else if(mixture_comp == 1) {
      true_sample[i, ] <- mvrnorm(n = 1, mu = 0 - mu, Sigma = diag(d))
    }
  }
  
  # Vector of subsampled test statistics
  ts_vec <- rep(NA, B)
  
  # Repeatedly subsample to get test statistic
  for(b in 1:B) {
    
    # Split Y into Y_0 and Y_1
    Y_0_indices <- sample(1:n_obs, size = floor(p * n_obs))
    
    Y_1_indices <- setdiff(1:n_obs, Y_0_indices)
    
    Y_0 <- matrix(true_sample[Y_0_indices, ], ncol = d)
    
    Y_1 <- matrix(true_sample[Y_1_indices, ], ncol = d)
    
    # Evaluate true density on D_0
    eval_true_D0 <- 
      (1 - p) * dmvnorm(x = Y_0, mean = rep(0, d), sigma = diag(d)) +
      p * dmvnorm(x = Y_0, mean = 0 - mu, sigma = diag(d))
    
    # Get log-concave MLE on D_0
    log_concave_D0 <- mlelcd(x = Y_0)
    
    # Evaluate log-concave MLE on D_0
    eval_log_concave_D0 <- dlcd(x = Y_0, lcd = log_concave_D0)
    
    # Compute test statistic
    ts_vec[b] <- exp(sum(log(eval_true_D0)) - sum(log(eval_log_concave_D0)))
    
    # Stop early if likelihood ratio already guarantees rejection
    if(compute_ts == 0 & sum(ts_vec[1:b]) > B / alpha) {
      ts_vec[(b + 1):B] <- 0
      break
    }
  }
  
  # Get average test statistic
  avg_ts <- mean(ts_vec)
  
  # Store average test statistic
  if(compute_ts == 1) {
    results[row, test_stat := avg_ts]
  }
  
  # Reject H_0 if test_stat >= 1/alpha
  results[row, reject := as.numeric(avg_ts >= 1/alpha)]
  
}

# Save simulation results
fwrite(results, file = paste0("sim_data/fig04_full_oracle_ddim_", line_number,
                              ".csv"))