# Benchmark permutation test
# Implement permutation test from Cule et al. (2010) to test
# H_0: log-concave versus H_1: not log-concave.
# Use two-component normal location model.
# Note: To load LogConcDEAD package, make sure XQuartz is running.

suppressMessages(library(LogConcDEAD))
suppressMessages(library(MASS))
suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(tidyverse))

# Read in arguments for file with all parameters and 
# line number for parameters for current simulation.

parameter_file <- "sim_params/tab01_perm_test_params.csv"
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
# equal_space_mu (indicator for mu structure. 
#                 If 1, each component is ||mu||*d^(-1/2).
#                 If 0, mu = (||mu||, 0, ... 0).
parameters <- parameter_df %>% slice(line_number)
d <- parameters$d
mu_norm <- parameters$mu_norm
n_obs <- parameters$n_obs
start_sim <- parameters$start_sim
n_sim <- parameters$n_sim
B <- parameters$B
equal_space_mu <- parameters$equal_space_mu

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
results <- data.table(Method = "Permutation test",
                      n_obs = n_obs, d = d, mu_norm = mu_norm, 
                      equal_space_mu = equal_space_mu, B = B,
                      sim = start_sim:(start_sim + n_sim - 1),
                      time_sec = NA_real_)

# Run simulations
for(row in 1:nrow(results)) {
  
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
  
  # Get start time for simulation
  start_time <- proc.time()[3]
  
  # Get MLE log-concave density
  dens_est <- mlelcd(x = true_sample)
  
  # Sample n_obs from MLE log-concave density
  dens_sample <- rlcd(n = n_obs, lcd = dens_est)
  
  # Combine all observations
  all_obs <- rbind(true_sample, dens_sample)

  # Initialize original sample test statistic
  orig_ts <- 0
  
  # Compute original sample test statistic
  for(i in 1:nrow(all_obs)) {
    
    # Get distance between all observations in true_sample and obs i
    P_n <- apply(matrix(true_sample, ncol = d), MARGIN = 1, 
                 FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))
    
    # Get distance between all observations in dens_sample and obs i
    P_n_star <- apply(matrix(dens_sample, ncol = d), MARGIN = 1, 
                      FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))
    
    # Get proportion of observations in sphere w/ each radius of size dist.
    # Save test stat that maximizes absolute difference.
    for(dist in unique(c(P_n, P_n_star))) {
      ts <- abs(mean(P_n <= dist) - mean(P_n_star <= dist))
      if(ts > orig_ts) {
        orig_ts <- ts
      }
    }
    
  }
  
  # Initialize permutation test statistics
  shuffle_ts <- rep(0, B)
  
  # Run permutation test
  for(b in 1:B) {
    
    # Shuffle order of observations
    index_shuffle <- sample(1:nrow(all_obs), size = nrow(all_obs), 
                            replace = FALSE)
    
    # Compute test statistic for permutation
    for(i in 1:nrow(all_obs)) {
      
      # Get distance between all observations in first sample and obs i
      P_n <- apply(matrix(all_obs[index_shuffle[1:n_obs], ], ncol = d), 
                   MARGIN = 1, 
                   FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))
      
      # Get distance between all observations in second sample and obs i
      P_n_star <- apply(matrix(all_obs[index_shuffle[(n_obs+1):(2*n_obs)], ], 
                               ncol = d),
                        MARGIN = 1, 
                        FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))
      
      # Get proportion of observations in sphere w/ each radius of size dist.
      # Save test stat that maximizes absolute difference.
      for(dist in unique(c(P_n, P_n_star))) {
        ts <- abs(mean(P_n <= dist) - mean(P_n_star <= dist))
        if(ts > shuffle_ts[b]) {
          shuffle_ts[b] <- ts
        }
      }
      
    }

  }
  
  # Get end time for simulation
  end_time <- proc.time()[3]
  
  # Store run time for simulation
  results[row, time_sec := end_time - start_time]

}

# Save simulation results
fwrite(results, file = paste0("sim_data/tab01_perm_test_", line_number, ".csv"))