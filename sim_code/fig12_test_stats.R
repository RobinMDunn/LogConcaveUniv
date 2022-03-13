# Get shuffled test statistics for permutation test from Cule et al. (2010).
# Testing H_0: log-concave versus H_1: not log-concave.
# Use two-component normal location model.
# Note: To load LogConcDEAD package, make sure XQuartz is running.

# Read in library
library(LogConcaveUniv)

# Read in arguments for file with all parameters and
# line number for parameters for current simulation.

parameter_file <- "sim_params/fig12_test_stats_params.csv"
line_number <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  parameter_file <- args[1]
  line_number <- as.numeric(args[2])
}

parameter_df <- data.table::fread(parameter_file)

# Assign arguments based on input.
# Arguments are d (dimension), mu_norm, n_obs (number of obs),
# sim_id (ID to keep track of simulations corresponding to single sample),
# and B (number of subsamples).

parameters <- parameter_df %>% dplyr::slice(line_number)
d <- parameters$d
mu_norm <- parameters$mu_norm
n_obs <- parameters$n_obs
sim_id <- parameters$sim_id
B <- parameters$B

# Create mu vector: (mu_norm, 0, ..., 0)
mu <- rep(0, d)
mu[1] <- mu_norm

# Proportions for the two Gaussian components
p <- 0.5

# Alpha level
alpha <- 0.1

# Create data frame to store results
results <- data.table::data.table(n_obs = n_obs, d = d, mu_norm = mu_norm,
                                  B = B, sim_id = sim_id,
                                  orig_test_stat = NA_real_,
                                  shuffle_test_stat = rep(NA_real_, B),
                                  reject = NA_real_)

# Generate sample from two-component normal location model
true_sample <- matrix(NA, nrow = n_obs, ncol = d)

for(i in 1:n_obs) {
  mixture_comp <- rbinom(n = 1, size = 1, prob = p)
  if(mixture_comp == 0) {
    true_sample[i, ] <- rnorm(n = d, mean = 0, sd = 1)
  } else if(mixture_comp == 1) {
    true_sample[i, ] <- MASS::mvrnorm(n = 1, mu = 0 - mu, Sigma = diag(d))
  }
}

# Run permutation test
perm_test_out <- LogConcaveUniv::permutation_test(data = true_sample,
                                                  B = B, alpha = alpha)

# Store original sample test statistic
results[, orig_test_stat := perm_test_out$orig_ts]

# Get shuffle permutation test statistics
for(b in 1:B) {

  # Store test statistic for this shuffle
  results[b, shuffle_test_stat := perm_test_out$shuffle_ts[b]]

}

# Reject H_0 if orig_ts > (B+1)*(1-alpha) quantile of shuffle_ts
results[, reject := perm_test_out$reject_null]

# Save simulation results
data.table::fwrite(results, file = paste0("sim_data/fig12_perm_test_stats_",
                                          line_number, ".csv"))
