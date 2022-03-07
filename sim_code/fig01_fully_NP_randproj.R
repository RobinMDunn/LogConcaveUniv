# Test H_0: log-concave versus H_1: not log-concave using universal LRT.
# Choose a random point along surface of d-dimensional unit sphere. 
# Project observations onto this vector, and rotate to x-axis.
# Check whether this one-dimensional projection is log-concave.
# Split data into D_0 and D_1. 
# Get kernel density estimate on D_1. Get log-concave MLE on D_0.
# Evaluate at all values in D_0.
# Repeat at multiple projections. Reject if mean test stat >= 1/alpha.

suppressMessages(library(logcondens))
suppressMessages(library(MASS))
suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(tidyverse))
suppressMessages(library(ks))

# Read in arguments for file with all parameters and 
# line number for parameters for current simulation.

parameter_file <- "sim_params/fig01_fully_NP_randproj_params.csv"
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
# n_proj (number of random projections),
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
n_proj <- parameters$n_proj
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
                      equal_space_mu = equal_space_mu, n_proj = n_proj, B = B,
                      sim = start_sim:(start_sim + n_sim - 1), alpha = alpha,
                      p_0 = p, avg_ts = NA_real_, reject = NA_real_)

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results), 
                       clear = T, show_after = 0)

# Run simulations to check whether to reject H_0
for(row in 1:nrow(results)) {
  
  #Increment progression bar
  pb$tick()
  
  # Generate sample from two-component normal location model
  true_sample <- matrix(NA, nrow = n_obs, ncol = d)
  
  for(i in 1:n_obs) {
    mixture_comp <- rbinom(n = 1, size = 1, prob = p)
    if(mixture_comp == 0) {
      true_sample[i, ] <- mvrnorm(n = 1, mu = rep(0, d), Sigma = diag(d))
    } else if(mixture_comp == 1) {
      true_sample[i, ] <- mvrnorm(n = 1, mu = 0 - mu, Sigma = diag(d))
    }
  }
  
  # Matrix of subsampled test statistics
  ts_mat <- matrix(NA, nrow = n_proj, ncol = B)
  
  # Repeatedly get test statistic on different projections
  for(proj in 1:n_proj) {
    
    # Get random vector for projection
    random_vector <- rnorm(n = d, mean = 0, sd = 1)
    
    random_vector <- random_vector / sqrt(sum(random_vector^2))
    
    for(b in 1:B) {
      
      # Split Y into Y_0 and Y_1
      Y_0_indices <- sample(1:n_obs, size = n_obs/2)
      
      Y_1_indices <- setdiff(1:n_obs, Y_0_indices)
      
      Y_0 <- matrix(true_sample[Y_0_indices, ], ncol = d)
      
      Y_1 <- matrix(true_sample[Y_1_indices, ], ncol = d)
      
      # Get projections of Y_0 and Y_1
      Y_0 <- as.numeric(Y_0 %*% random_vector)
      
      Y_1 <- as.numeric(Y_1 %*% random_vector)
      
      # Get KDE estimate on D_1.
      kde_D1 <- ks::kde(x = Y_1, eval.points = Y_0,
                        h = hpi(x = Y_1), binned = FALSE)
      
      # Evaluate KDE on D_0
      eval_kde_D0 <- kde_D1$estimate
      
      # Get log-concave MLE on D_0
      log_concave_D0 <- logConDens(x = Y_0, smoothed = F)

      # Evaluate log-concave MLE on D_0
      eval_log_concave_D0 <- evaluateLogConDens(
        xs = Y_0, res = log_concave_D0, which = 2)[, 3]
      
      # Store test stat for the projection
      ts_mat[proj, b] <- 
        exp(sum(log(eval_kde_D0)) - sum(log(eval_log_concave_D0)))
      
      # Stop if test stat is NA
      stopifnot(!is.na(ts_mat[proj, b]))
      
    }
    
    # Stop early if likelihood ratio already guarantees rejection
    if(compute_ts == 0 & sum(ts_mat, na.rm = T) >= n_proj * B / alpha) {
      ts_mat[is.na(ts_mat)] <- 0
      break
    }

  }
  
  # Get average subsampled test stat across projections
  avg_test_stat <- mean(ts_mat)
  
  # Store average test stat
  if(compute_ts == 1) {
    results[row, avg_ts := avg_test_stat]
  }

  # Reject H_0 if avg_ts >= 1/alpha
  results[row, reject := as.numeric(avg_test_stat >= 1/alpha)]
 
}

# Save simulation results
fwrite(results, file = paste0("sim_data/fig01_fully_NP_randproj_",
                              line_number, ".csv"))
