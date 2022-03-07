# Implement permutation test from Cule et al. (2010) to test
# H_0: log-concave versus H_1: not log-concave.
# Use beta distribution. Density is log-concave iff alpha_param, beta_param >= 1.
# Note: Use logcondens since d = 1.

suppressMessages(library(logcondens))
suppressMessages(library(MASS))
suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(tidyverse))

# Read in arguments for file with all parameters and 
# line number for parameters for current simulation.

parameter_file <- "sim_params/fig17_perm_test_params.csv"
line_number <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  parameter_file <- args[1]
  line_number <- as.numeric(args[2])
}

parameter_df <- fread(parameter_file)

# Assign arguments based on input.
# Arguments are alpha_param (Beta dist param 1), beta_param (Beta dist param 2), 
# n_obs (number of obs),
# start_sim (index of starting sim), n_sim (number of sims), 
# B (number of subsamples)
parameters <- parameter_df %>% slice(line_number)
alpha_param <- parameters$alpha_param
beta_param <- parameters$beta_param
n_obs <- parameters$n_obs
start_sim <- parameters$start_sim
n_sim <- parameters$n_sim
B <- parameters$B

# Alpha level
alpha_level <- 0.1

# Create data frame to store results
results <- data.table(n_obs = n_obs, d = 1, 
                      sim = start_sim:(start_sim + n_sim - 1),
                      alpha_param = alpha_param,
                      beta_param = beta_param, B = B,
                      reject = NA_real_)

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results) * B, 
                       clear = T, show_after = 0)

# Run simulations
for(row in 1:nrow(results)) {
  
  # Extract shape1 (alpha) and shape2 (beta) params
  alpha_param <- results[row, alpha_param]
  beta_param <- results[row, beta_param]
  
  # Generate sample from beta density
  true_sample <-
    matrix(rbeta(n = n_obs, shape1 = alpha_param, shape2 = beta_param), ncol = 1)
  
  # Get MLE log-concave density. Use logcondens since d = 1.
  dens_est <- logConDens(x = true_sample, smoothed = F)
  
  # Sample n_obs from MLE log-concave density
  dens_sample <- rlogcon(n = n_obs, x0 = sort(true_sample))$X
  
  # Combine all observations
  all_obs <- rbind(true_sample, as.matrix(dens_sample, ncol = 1))

  # Initialize original sample test statistic
  orig_ts <- 0
  
  # Compute original sample test statistic
  for(i in 1:nrow(all_obs)) {
    
    # Get distance between all observations in true_sample and obs i
    P_n <- apply(matrix(true_sample, ncol = 1), MARGIN = 1, 
                 FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))
    
    # Get distance between all observations in dens_sample and obs i
    P_n_star <- apply(matrix(dens_sample, ncol = 1), MARGIN = 1, 
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
    
    # Update progress bar
    pb$tick()
    
    # Shuffle order of observations
    index_shuffle <- sample(1:nrow(all_obs), size = nrow(all_obs), 
                            replace = FALSE)
    
    # Compute test statistic for permutation
    for(i in 1:nrow(all_obs)) {
      
      # Get distance between all observations in first sample and obs i
      P_n <- apply(matrix(all_obs[index_shuffle[1:n_obs], ], ncol = 1), 
                   MARGIN = 1, 
                   FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))
      
      # Get distance between all observations in second sample and obs i
      P_n_star <- apply(matrix(all_obs[index_shuffle[(n_obs+1):(2*n_obs)], ], 
                               ncol = 1),
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
  
  # Reject H_0 if orig_ts > (B+1)*(1-alpha) quantile of shuffle_ts
  results[row, reject := as.numeric(orig_ts > sort(shuffle_ts)[ceiling((B+1)*(1-alpha_level))])]

}

# Save simulation results
fwrite(results, file = paste0("sim_data/fig17_perm_test_", line_number, ".csv"))