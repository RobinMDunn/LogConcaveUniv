# Implement permutation test from Cule et al. (2010) to test
# H_0: log-concave versus H_1: not log-concave.
# Use two-component normal location model.
# Note: To load LogConcDEAD package, make sure XQuartz is running.

# Read in library
library(LogConcaveUniv)

# Read in arguments for file with all parameters and
# line number for parameters for current simulation.

parameter_file <- "sim_params/fig03_perm_test_params.csv"
line_number <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  parameter_file <- args[1]
  line_number <- as.numeric(args[2])
}

parameter_df <- data.table::fread(parameter_file)

# Assign arguments based on input.
# Arguments are d (dimension), mu_norm, n_obs (number of obs),
# start_sim (index of starting sim), n_sim (number of sims),
# B (number of subsamples),
# equal_space_mu (indicator for mu structure in second component.
#                 If 1, each component is ||mu||*d^(-1/2).
#                 If 0, mu = (||mu||, 0, ... 0).
parameters <- parameter_df %>% dplyr::slice(line_number)
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
results <- data.table::data.table(n_obs = n_obs, d = d, mu_norm = mu_norm,
                                  equal_space_mu = equal_space_mu, B = B,
                                  sim = start_sim:(start_sim + n_sim - 1),
                                  reject = NA_real_)

# Set up progress bar
pb <- progress::progress_bar$new(format = "sim :current / :total [:bar] :eta",
                                 total = nrow(results),
                                 clear = T, show_after = 0)

# Run simulations
for(row in 1:nrow(results)) {

  # Update progress bar
  pb$tick()

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

  # Run permutation test to determine whether to reject H_0
  reject_null <- LogConcaveUniv::permutation_test(data = true_sample, B = B,
                                                  alpha = alpha)$reject_null

  results[row, reject := reject_null]

}

# Save simulation results
data.table::fwrite(results, file = paste0("sim_data/fig03_perm_test_",
                                          line_number, ".csv"))
