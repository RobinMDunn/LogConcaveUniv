# Test H_0: log-concave versus H_1: not log-concave using universal LRT.
# Split data into D_0 and D_1.
# Use beta distribution. Density is log-concave iff alpha_param, beta_param >= 1.
# Get bounded KDE estimate on D_1, using plug-in bw. Get log-concave MLE on D_0.
# Evaluate at all values in D_0.

# Read in library
library(LogConcaveUniv)

# Read in arguments for file with all parameters and
# line number for parameters for current simulation.

parameter_file <- "sim_params/fig17_fully_NP_params.csv"
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
results <-
  data.table::data.table(n_obs = n_obs, d = 1,
                         sim = start_sim:(start_sim + n_sim - 1),
                         alpha_param = alpha_param,
                         beta_param = beta_param, B = B,
                         test_stat = NA_real_,
                         reject = NA_real_)

# Set up progress bar
pb <- progress::progress_bar$new(format = "sim :current / :total [:bar] :eta",
                                 total = nrow(results) * B,
                                 clear = T, show_after = 0)

# Run simulations to check whether to reject H_0
for(row in 1:nrow(results)) {

  # Extract shape1 (alpha) and shape2 (beta) params
  alpha_param <- results[row, alpha_param]
  beta_param <- results[row, beta_param]

  # Generate sample from beta density
  true_sample <-
    matrix(rbeta(n = n_obs, shape1 = alpha_param, shape2 = beta_param), ncol = 1)

  # Vector of subsampled test statistics
  ts_vec <- rep(NA, B)

  # Repeatedly subsample to get test statistic
  for(b in 1:B) {

    #Increment progression bar
    pb$tick()

    # Split Y into Y_0 and Y_1
    Y_0_indices <- sample(1:n_obs, size = n_obs / 2)

    Y_1_indices <- setdiff(1:n_obs, Y_0_indices)

    Y_0 <- matrix(true_sample[Y_0_indices, ], ncol = 1)

    Y_1 <- matrix(true_sample[Y_1_indices, ], ncol = 1)

    # Get KDE on D_1
    kde_D1 <- kde1d::kde1d(x = Y_1, xmin = 0, xmax = 1)

    # Evaluate KDE on D_0
    eval_kde_D0 <- kde1d::dkde1d(x = Y_0, obj = kde_D1)

    # Get log-concave MLE on D_0
    log_concave_D0 <- logcondens::logConDens(x = Y_0, smoothed = FALSE)

    # Evaluate log-concave MLE on D_0
    eval_log_concave_D0 <-
      logcondens::evaluateLogConDens(xs = Y_0, res = log_concave_D0)[, 3]

    # Compute test statistic
    ts_vec[b] <- exp(sum(log(eval_kde_D0)) - sum(log(eval_log_concave_D0)))
  }

  # Get average test statistic
  avg_ts <- mean(ts_vec)

  # Store average test statistic
  results[row, test_stat := avg_ts]

  # Reject H_0 if test_stat >= 1/alpha
  results[row, reject := as.numeric(avg_ts >= 1/alpha_level)]

}

# Save simulation results
data.table::fwrite(results, file = paste0("sim_data/fig17_fully_NP_",
                                          line_number, ".csv"))
