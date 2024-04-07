# Get densities to plot appearance of log-concave MLEs from test of
# H_0: log-concave versus H_1: not log-concave using universal LRT.
# About 1 hour to run mlelcd log-concave MLE function on d = 1, n = 5000.
# Less than 5 seconds to run logcondens on d = 1, n = 5000.

# Read in library
library(LogConcaveUniv)

# Proportions for the two Gaussian components
p <- 0.5

# Create data frame to store density_df
density_df <-
  data.table::data.table(d = 1,
                         expand.grid(x = seq(-7, 3, length.out = 201),
                                     y = NA_real_,
                                     mu_norm = c(0, 2, 4)),
                         n_obs = 5000,
                         true_density = NA_real_,
                         LogConcDEAD_density = NA_real_,
                         logcondens_density = NA_real_,
                         mean_loglik_true_dens = NA_real_,
                         mean_loglik_LogConcDEAD = NA_real_,
                         mean_loglik_logcondens = NA_real_)

# Set simulation numbers for seed values
density_df[mu_norm == 0 & n_obs == 5000 & d == 1, sim := 10]
density_df[mu_norm == 2 & n_obs == 5000 & d == 1, sim := 11]
density_df[mu_norm == 4 & n_obs == 5000 & d == 1, sim := 12]

# Set up progress bar
pb <- progress::progress_bar$new(format = "sim :current / :total [:bar] :eta",
                                 total = length(unique(density_df$sim)),
                                 clear = T, show_after = 0)

# Run simulations to check whether to reject H_0
for(sim_val in 10:12) {

  # Update progress bar
  pb$tick()

  # Set seed
  set.seed(sim_val)

  # Extract dimension
  d <- density_df[sim == sim_val, d][1]

  # Check that d = 1
  stopifnot(d == 1)

  # Create mu vector: (mu_norm, 0, ..., 0)
  mu_norm <- density_df[sim == sim_val, mu_norm][1]
  mu <- rep(0, d)
  mu[1] <- mu_norm

  # Extract n_obs
  n_obs <- density_df[sim == sim_val, n_obs][1]

  # Generate sample from two-component normal location model
  true_sample <- rep(NA, length = n_obs)

  for(i in 1:n_obs) {
    mixture_comp <- rbinom(n = 1, size = 1, prob = p)
    if(mixture_comp == 0) {
      true_sample[i] <- rnorm(n = d, mean = 0, sd = 1)
    } else if(mixture_comp == 1) {
      true_sample[i] <- rnorm(n = 1, mean = 0 - mu, sd = 1)
    }
  }

  # Get MLE log-concave estimate using LogConcDEAD package
  mle_LogConcDEAD <- LogConcDEAD::mlelcd(x = true_sample)

  # Get MLE log-concave estimate using logcondens package
  mle_logcondens <- logcondens::logConDens(x = true_sample)

  # Extract x and y coords
  x_coord <- density_df[sim == sim_val, x]

  xy_mat <- matrix(x_coord, ncol = 1)

  # True values for Gaussian mixture
  true_mu <- matrix(rbind(0 - mu, rep(0, d)), nrow = 2)
  true_sigma <- array(dim = c(d, d, 2))
  true_sigma[, , 1] <- true_sigma[, , 2] <- diag(d)

  # Evaluate true Gaussian mixture density
  density_df[sim == sim_val, true_density := EMCluster::dmixmvn(x = xy_mat,
                            pi = c(0.5, 0.5), Mu = true_mu,
                            LTSigma = EMCluster::variance2LTSigma(true_sigma))]

  # Evaluate MLE log-concave density using LogConcDEAD package
  density_df[sim == sim_val,
             LogConcDEAD_density := LogConcDEAD::dlcd(x = xy_mat,
                                                      lcd = mle_LogConcDEAD)]

  # Evaluate MLE log-concave estimate using logcondens package
  density_df[sim == sim_val,
             logcondens_density := logcondens::evaluateLogConDens(
               xs = xy_mat, res = mle_logcondens, which = 2)[, 3]]

  # Get (1/n)*loglik on observed data, true Gaussian mixture density
  mean_loglik_true <-
    mean(EMCluster::dmixmvn(x = matrix(true_sample, ncol = 1), pi = c(0.5, 0.5),
                 Mu = true_mu,
                 LTSigma = EMCluster::variance2LTSigma(true_sigma), log = TRUE))

  density_df[sim == sim_val, mean_loglik_true_dens := mean_loglik_true]

  # Get (1/n)*loglik on observed data, LogConcDEAD density
  mean_loglik_LogConcDEAD_eval <-
    mean(LogConcDEAD::dlcd(x = true_sample, lcd = mle_LogConcDEAD,
                           uselog = TRUE))

  density_df[sim == sim_val,
             mean_loglik_LogConcDEAD := mean_loglik_LogConcDEAD_eval]

  # Get (1/n)*loglik on observed data, logcondens density
  mean_loglik_logcondens_eval <- mean(logcondens::evaluateLogConDens(
    xs = true_sample, res = mle_logcondens, which = 1)[, 2])

  density_df[sim == sim_val,
             mean_loglik_logcondens := mean_loglik_logcondens_eval]

}

# Save simulation density_df
data.table::fwrite(density_df, file = "sim_data/fig02_densities.csv")
