# Benchmark partial oracle, d-dim
# Test H_0: log-concave versus H_1: not log-concave using universal LRT.
# Split data into D_0 and D_1. 
# Get KDE estimate on D_1. Get log-concave MLE on D_0.
# Evaluate at all values in D_0.

suppressMessages(library(LogConcDEAD))
suppressMessages(library(MASS))
suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(tidyverse))
suppressMessages(library(mclust))

# Read in arguments for file with all parameters and 
# line number for parameters for current simulation.

parameter_file <- "sim_params/tab01_partial_oracle_ddim_params.csv"
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
results <- data.table(Method = "Partial oracle, d-dim",
                      n_obs = n_obs, d = d, mu_norm = mu_norm, 
                      equal_space_mu = equal_space_mu, B = B,
                      sim = start_sim:(start_sim + n_sim - 1), alpha = alpha,
                      p_0 = p, time_sec = NA_real_)

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
  
  # Get start time for simulation
  start_time <- proc.time()[3]
  
  # Vector of subsampled test statistics
  ts_vec <- rep(NA, B)
  
  # Repeatedly subsample to get test statistic
  for(b in 1:B) {
    
    # Split Y into Y_0 and Y_1
    Y_0_indices <- sample(1:n_obs, size = floor(p * n_obs))
    
    Y_1_indices <- setdiff(1:n_obs, Y_0_indices)
    
    Y_0 <- matrix(true_sample[Y_0_indices, ], ncol = d)
    
    Y_1 <- matrix(true_sample[Y_1_indices, ], ncol = d)
    
    # Remove previously fitted mclust_dens_D1 object
    if(exists("mclust_dens_D1")) { rm(mclust_dens_D1) }
    
    # Get two-component Gaussian mixture estimate on D_1.
    # modelNames:
    #   - V: unequal variance
    #   - E: equal variance
    #   - VVV: unstructured (other than symmetric) covariance matrices.
    #   - VEV: varying volume, equal shape, varying orientation (Banfield & Raftery, 1993)
    if(d == 1) {
     
      mclust_dens_D1 <- try(densityMclust(data = Y_1, G = 2, modelNames = "V",
                                          warn = FALSE, verbose = FALSE))
      if(is(mclust_dens_D1, "try-error") | is.null(mclust_dens_D1)) {
        mclust_dens_D1 <- try(densityMclust(data = Y_1, G = 2, modelNames = "E",
                                            warn = FALSE, verbose = FALSE))
      }
      
    } else {
      
      mclust_dens_D1 <- try(densityMclust(data = Y_1, G = 2, modelNames = "VVV", 
                                          warn = FALSE, verbose = FALSE))
      if(is(mclust_dens_D1, "try-error") | is.null(mclust_dens_D1)) {
        mclust_dens_D1 <- try(densityMclust(data = Y_1, G = 2, modelNames = "VEV", 
                                            warn = FALSE, verbose = FALSE))
      }
      
    }
    
    # If mclust_dens_D1 is not an Mclust object, fit a single Normal density.
    # (This is very rare in simulations.)
    if(!("Mclust" %in% class(mclust_dens_D1))) {
      print("Fitting single Normal density")
      if(d == 1) {
        mclust_dens_D1 <- try(densityMclust(data = Y_1, G = 1, modelNames = "V",
                                            warn = FALSE, verbose = FALSE,
                                            plot = FALSE))
      } else {
        mclust_dens_D1 <- try(densityMclust(data = Y_1, G = 1, modelNames = "VVV", 
                                            warn = FALSE, verbose = FALSE,
                                            plot = FALSE))
      }
    }
    
    # # Parameters of Gaussian mixture
    # prop_est_D1 <- mclust_dens_D1$parameters$pro
    # mean_est_D1 <- mclust_dens_D1$parameters$mean (for 2 classes)
    # var_est_D1 <- mclust_dens_D1$parameters$variance$sigma # array of 2 d x d matrices for d > 1
    # var_est_D1 <- mclust_dens_D1$parameters$variance$sigmasq # vector length 2 for d = 1
    
    # Evaluate Gaussian mixture on D_0
    eval_gauss_mix_D0 <- predict.densityMclust(object = mclust_dens_D1, 
                                               newdata = Y_0, what = "dens")
    
    # Stop if max evaluated density is over 100 (probably convergence issue)
    # (Never occurs in simulations.)
    stopifnot(max(eval_gauss_mix_D0) < 100)
    
    ## Same as
    # prop_est_D1[1]*mclust::dmvnorm(data = Y_0, mean = mean_est_D1[,1],
    #                                sigma = var_est_D1[,,1]) +
    #   prop_est_D1[2]*mclust::dmvnorm(data = Y_0, mean = mean_est_D1[,2],
    #                                  sigma = var_est_D1[,,2])
    
    # Get log-concave MLE on D_0
    log_concave_D0 <- mlelcd(x = Y_0)
    
    # Evaluate log-concave MLE on D_0
    eval_log_concave_D0 <- dlcd(x = Y_0, lcd = log_concave_D0)
    
    # Compute test statistic
    ts_vec[b] <- exp(sum(log(eval_gauss_mix_D0)) - sum(log(eval_log_concave_D0)))
    
    # Stop early if likelihood ratio already guarantees rejection
    if(compute_ts == 0 & sum(ts_vec[1:b]) > B / alpha) {
      ts_vec[(b + 1):B] <- 0
      break
    }
  }
  
  # Get end time for simulation
  end_time <- proc.time()[3]
  
  # Store run time for simulation
  results[row, time_sec := end_time - start_time]
  
}

# Save simulation results
fwrite(results, file = paste0("sim_data/tab01_partial_oracle_ddim_",
                              line_number, ".csv"))
