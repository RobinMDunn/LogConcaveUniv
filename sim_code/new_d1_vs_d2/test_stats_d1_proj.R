# Test H_0: log-concave versus H_1: not log-concave using universal LRT.
# Choose a random point along surface of d-dimensional unit sphere.
# Project observations onto this vector, and rotate to x-axis.
# Check whether this one-dimensional projection is log-concave.
# Split data into D_0 and D_1.
# Get Gaussian mixture estimate on D_1. Get log-concave MLE on D_0.
# Evaluate at all values in D_0.
# True density is (1/2)N((0,0), I_2) + (1/2)N((-6,0), I_2).
# Store test statistics from partial oracle projection approach.

# Read in library
library(LogConcaveUniv)

# Set parameter values
d <- 2
mu_1 <- 6
mu_2 <- 0
n_obs <- 100
n_sim <- 1000
B <- 100
equal_space_mu <- 0
compute_ts <- 1
p <- 0.5
alpha <- 0.1

# Use a single projection per dataset for illustration
n_proj <- 1

# Choose number of cores over which to run simulation
# (may want to change based on computing infrastructure)
n_cores <- 50

# Create data frame to store results
results <- data.table::data.table(n_obs = n_obs, d = d,
                                  mu_1 = mu_1, mu_2 = mu_2,
                                  equal_space_mu = equal_space_mu,
                                  compute_ts = compute_ts, B = B,
                                  n_proj = n_proj, sim = 1:n_sim,
                                  alpha = alpha, p_0 = p)

# Function to run simulation
run_simulation_partial_oracle_randproj <- 
  function(n_obs, d, mu_1, mu_2, equal_space_mu, compute_ts, B, n_proj, sim, 
           alpha, p_0) {

  # Get random vector for projection
  # (typically done inside partial_oracle_randproj)
  random_vector <- rnorm(n = d, mean = 0, sd = 1)

  random_vector <- random_vector / sqrt(sum(random_vector^2))

  # Generate sample from two-component normal location model
  true_sample <- matrix(NA, nrow = n_obs, ncol = d)

  for(i in 1:n_obs) {
    mixture_comp <- rbinom(n = 1, size = 1, prob = p_0)
    if(mixture_comp == 0) {
        true_sample[i, ] <- rnorm(n = d, mean = 0, sd = 1)
    } else if(mixture_comp == 1) {
        true_sample[i, ] <- MASS::mvrnorm(n = 1, mu = 0 - c(mu_1, mu_2), 
                                          Sigma = diag(d))
    }
  }

  # Project data onto random_vector
  true_sample <- true_sample %*% random_vector

  # Run partial oracle random projection test to determine whether to reject H_0
  test_out <- 
    LogConcaveUniv::partial_oracle_randproj(data = true_sample, B = B,
                                            n_proj = n_proj, alpha = alpha,
                                            compute_ts = compute_ts)

  return(list(random_vector_1 = random_vector[1],
              random_vector_2 = random_vector[2],
              test_stat = test_out$test_stat,
              reject_null = test_out$reject_null))
}

# Run simulations and compute test statistics
sim_out <- clustermq::Q_rows(df = results,
                             fun = run_simulation_partial_oracle_randproj,
                             n_jobs = n_cores)

# Store results
results$random_vector_1 <- unlist(lapply(sim_out, FUN = function(x) x$random_vector_1))
results$random_vector_2 <- unlist(lapply(sim_out, FUN = function(x) x$random_vector_2))
results$test_stat <- unlist(lapply(sim_out, FUN = function(x) x$test_stat))
results$reject_null <- unlist(lapply(sim_out, FUN = function(x) x$reject_null))

# Save reuslts
data.table::fwrite(results,
                   file = "sim_data/new_d1_vs_d2/test_stats_d1_proj.csv")
