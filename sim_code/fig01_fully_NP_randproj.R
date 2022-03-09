# Test H_0: log-concave versus H_1: not log-concave using universal LRT.
# Choose a random point along surface of d-dimensional unit sphere.
# Project observations onto this vector, and rotate to x-axis.
# Check whether this one-dimensional projection is log-concave.
# Split data into D_0 and D_1.
# Get kernel density estimate on D_1. Get log-concave MLE on D_0.
# Evaluate at all values in D_0.
# Repeat at multiple projections. Reject if mean test stat >= 1/alpha.

# Read in library
library(LogConcaveUniv)

# Read in arguments for file with all parameters and
# line number for parameters for current simulation.
parameter_file <- "sim_params/fig01_fully_NP_randproj_params.csv"
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
# n_proj (number of random projections),
# B (number of subsamples),
# equal_space_mu (indicator for mu structure in second component.
#                 If 1, each component is ||mu||*d^(-1/2).
#                 If 0, mu = (||mu||, 0, ... 0),
# compute_ts (indicator for whether to compute test stat.
#             If 1, computes test stat.
#             If 0, does not compute test stat, allowing for early rejection.)

parameters <- parameter_df %>% dplyr::slice(line_number)
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
results <- data.table::data.table(
  n_obs = n_obs, d = d, mu_norm = mu_norm,
  equal_space_mu = equal_space_mu, n_proj = n_proj, B = B,
  sim = start_sim:(start_sim + n_sim - 1), alpha = alpha,
  p_0 = p, avg_ts = NA_real_, reject = NA_real_
)

# Set up progress bar
pb <- progress::progress_bar$new(format = "sim :current / :total [:bar] :eta",
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
      true_sample[i, ] <- MASS::mvrnorm(n = 1, mu = rep(0, d), Sigma = diag(d))
    } else if(mixture_comp == 1) {
      true_sample[i, ] <- MASS::mvrnorm(n = 1, mu = 0 - mu, Sigma = diag(d))
    }
  }

  # Run fully nonparametric random projection test to determine
  # whether to reject H_0
  test_out <- fully_NP_randproj(data = true_sample, B = B, n_proj = n_proj,
                                alpha = alpha, compute_ts = compute_ts)

  results[row, avg_ts := test_out$test_stat]
  results[row, reject := test_out$reject_null]

}

# Save simulation results
data.table::fwrite(results,
                   file = paste0("sim_data/fig01_fully_NP_randproj_",
                                 line_number, ".csv"))
