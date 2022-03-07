# Test H_0: log-concave versus H_1: not log-concave using universal LRT.
# Split data into D_0 and D_1. 
# Consider each component separately (d >= 2).
# Get kernel density estimate on D_1. Get log-concave MLE on D_0.
# Evaluate at all values in D_0.
# Reject if T_{n,k} >= d/alpha for at least one k \in {1, 2, ..., d}.

suppressMessages(library(logcondens))
suppressMessages(library(MASS))
suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(tidyverse))
suppressMessages(library(ks))

# Read in arguments for file with all parameters and 
# line number for parameters for current simulation.

parameter_file <- "sim_params/fig04_fully_NP_axis_params.csv"
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
results <- data.table(n_obs = n_obs, d = d, mu_norm = mu_norm,
                      equal_space_mu = equal_space_mu, B = B,
                      sim = start_sim:(start_sim + n_sim - 1), alpha = alpha,
                      p_0 = p, reject = NA_real_)

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results) * d, 
                       clear = T, show_after = 0)

# Run simulations to check whether to reject H_0
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
  
  # Matrix of subsampled test statistics
  ts_mat <- matrix(NA, nrow = B, ncol = d)
  
  # Consider each dimension separately
  for(d_val in 1:d) {
    
    #Increment progress bar
    pb$tick()
    
    # Repeatedly subsample to get subsampled test statistic on each dimension
    for(b in 1:B) {
      
      # Split Y into Y_0 and Y_1
      Y_0_indices <- sample(1:n_obs, size = n_obs/2)
      
      Y_1_indices <- setdiff(1:n_obs, Y_0_indices)
      
      Y_0 <- matrix(true_sample[Y_0_indices, ], ncol = d)
      
      Y_1 <- matrix(true_sample[Y_1_indices, ], ncol = d)
      
      # Get KDE estimate for each dimension on D_1.
      # kde_D1 <- ks::kde(x = Y_1[, d_val], eval.points = Y_0[, d_val])
      kde_D1 <- ks::kde(x = Y_1[, d_val], eval.points = Y_0[, d_val],
                        h = hpi(x = Y_1[, d_val]), binned = FALSE)
      
      # Evaluate KDE on D_0
      eval_kde_D0 <- kde_D1$estimate
      
      # Get log-concave MLE on D_0
      log_concave_D0 <- logConDens(x = Y_0[, d_val])
      
      # Evaluate log-concave MLE on D_0
      eval_log_concave_D0 <- evaluateLogConDens(
        xs = Y_0[, d_val], res = log_concave_D0, which = 2)[, 3]
      
      # Store dimension d, subsample b test stat
      ts_mat[b, d_val] <- 
        exp(sum(log(eval_kde_D0)) - sum(log(eval_log_concave_D0)))
      
      # Break if you would reject based on current info
      if(sum(apply(ts_mat, MARGIN = 2, FUN = function(x) sum(x, na.rm = T)) >= 
          B * d /alpha) >= 1 & b < B) {
        ts_mat[(b+1):B, d_val] <- 0
        break
      }
      
    }
    
  }
  
  # Get average subsampled test stat for each dimension
  test_stats <- apply(ts_mat, MARGIN = 2, FUN = mean)
  
  # Reject H_0 if some dimension's test_stat >= d/alpha
  results[row, 
          reject := as.numeric(sum(test_stats >= d/alpha) >= 1)]
 
}

# Save simulation results
fwrite(results, file = paste0("sim_data/fig04_fully_NP_axis_", line_number, ".csv"))
