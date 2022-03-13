# Implement permutation test from Cule et al. (2010) to test
# H_0: log-concave versus H_1: not log-concave.
# Use beta distribution. Density is log-concave iff alpha_param, beta_param >= 1.
# Note: Use logcondens since d = 1.

# Read in library
library(LogConcaveUniv)

# Read in arguments for file with all parameters and
# line number for parameters for current simulation.

parameter_file <- "sim_params/fig17_perm_test_params.csv"
line_number <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  parameter_file <- args[1]
  line_number <- as.numeric(args[2])
}

parameter_df <- data.table::fread(parameter_file)

# Assign arguments based on input.
# Arguments are alpha_param (Beta dist param 1), beta_param (Beta dist param 2),
# n_obs (number of obs),
# start_sim (index of starting sim), n_sim (number of sims),
# B (number of subsamples)
parameters <- parameter_df %>% dplyr::slice(line_number)
alpha_param <- parameters$alpha_param
beta_param <- parameters$beta_param
n_obs <- parameters$n_obs
start_sim <- parameters$start_sim
n_sim <- parameters$n_sim
B <- parameters$B

# Alpha level
alpha_level <- 0.1

# Create data frame to store results
results <- data.table::data.table(n_obs = n_obs, d = 1,
                                  sim = start_sim:(start_sim + n_sim - 1),
                                  alpha_param = alpha_param,
                                  beta_param = beta_param, B = B,
                                  reject = NA_real_)

# Set up progress bar
pb <- progress::progress_bar$new(format = "sim :current / :total [:bar] :eta",
                                 total = nrow(results),
                                 clear = T, show_after = 0)

# Run simulations
for(row in 1:nrow(results)) {

  # Update progress bar
  pb$tick()

  # Extract shape1 (alpha) and shape2 (beta) params
  alpha_param <- results[row, alpha_param]
  beta_param <- results[row, beta_param]

  # Generate sample from beta density
  true_sample <-
    matrix(rbeta(n = n_obs, shape1 = alpha_param, shape2 = beta_param), ncol = 1)

  # Run permutation test to determine whether to reject H_0
  reject_null <- LogConcaveUniv::permutation_test(data = true_sample, B = B,
                                                  alpha = alpha)$reject_null

  results[row, reject := reject_null]

}

# Save simulation results
data.table::fwrite(results, file = paste0("sim_data/fig17_perm_test_",
                                          line_number, ".csv"))
