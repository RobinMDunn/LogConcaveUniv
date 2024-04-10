# Get densities to plot appearance of log-concave MLEs from test of
# H_0: log-concave versus H_1: not log-concave using universal LRT.
# Less than 1 second to run mlelcd log-concave MLE function on d = 2, n = 50.

# Read in library
library(LogConcaveUniv)

# Proportions for the two Gaussian components
p <- 0.5

# Create data frame to store density_df
density_df <-
  data.table::data.table(d = 2,
                         expand.grid(x = seq(-6, 3, length.out = 201),
                                     y = seq(-6, 3, length.out = 201),
                                     mu_norm = c(0, 2, 4)),
                         n_obs = 50,
                         true_density = NA_real_,
                         LogConcDEAD_density = NA_real_,
                         logcondens_density = NA_real_,
                         mean_loglik_true_dens = NA_real_,
                         mean_loglik_LogConcDEAD = NA_real_,
                         mean_loglik_logcondens = NA_real_)

# Set simulation numbers for seed values
density_df[mu_norm == 0 & n_obs == 50 & d == 2, sim := 1]
density_df[mu_norm == 2 & n_obs == 50 & d == 2, sim := 2]
density_df[mu_norm == 4 & n_obs == 50 & d == 2, sim := 3]

# Data frame to store simulated points
density_points <-
  data.table::data.table(mu_norm = c(rep(c("||u|| = 0", "||u|| = 2",
                                           "||u|| = 4"), each = 50)),
                         sim = rep(c(1, 2, 3), each = 50),
                         x = NA_real_,
                         y = NA_real_)

# Set up progress bar
pb <- progress::progress_bar$new(format = "sim :current / :total [:bar] :eta",
                                 total = length(unique(density_df$sim)),
                                 clear = T, show_after = 0)

# Run simulations to check whether to reject H_0
for(sim_val in 1:3) {

  # Update progress bar
  pb$tick()

  # Set seed
  set.seed(sim_val)

  # Extract dimension
  d <- density_df[sim == sim_val, d][1]

  # Create mu vector: (mu_norm, 0, ..., 0)
  mu_norm <- density_df[sim == sim_val, mu_norm][1]
  mu <- rep(0, d)
  mu[1] <- mu_norm

  # Extract n_obs
  n_obs <- density_df[sim == sim_val, n_obs][1]

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

  # Store simulated points
  density_points[sim == sim_val, x := true_sample[, 1]]
  density_points[sim == sim_val, y := true_sample[, 2]]

  # Get MLE log-concave estimate using LogConcDEAD package
  mle_LogConcDEAD <- LogConcDEAD::mlelcd(x = true_sample)

  # Extract x and y coords
  x_coord <- density_df[sim == sim_val, x]

  y_coord <- density_df[sim == sim_val, y]
  xy_mat <- cbind(x_coord, y_coord)

  # True values for Gaussian mixture
  true_mu <- matrix(rbind(0 - mu, rep(0, d)), nrow = 2)
  true_sigma <- array(dim = c(d, d, 2))
  true_sigma[, , 1] <- true_sigma[, , 2] <- diag(d)

  # Evaluate true Gaussian mixture density
  density_df[sim == sim_val,
             true_density := EMCluster::dmixmvn(
               x = xy_mat, pi = c(0.5, 0.5), Mu = true_mu,
               LTSigma = EMCluster::variance2LTSigma(true_sigma))]

  # Evaluate MLE log-concave density using LogConcDEAD package
  density_df[sim == sim_val,
             LogConcDEAD_density := LogConcDEAD::dlcd(x = xy_mat,
                                                      lcd = mle_LogConcDEAD)]

  # Get (1/n)*loglik on observed data, true Gaussian mixture density
  mean_loglik_true <-
    mean(EMCluster::dmixmvn(
      x = true_sample, pi = c(0.5, 0.5), Mu = true_mu,
      LTSigma = EMCluster::variance2LTSigma(true_sigma), log = TRUE))

  density_df[sim == sim_val, mean_loglik_true_dens := mean_loglik_true]

  # Get (1/n)*loglik on observed data, LogConcDEAD density
  mean_loglik_LogConcDEAD_eval <-
    mean(LogConcDEAD::dlcd(x = true_sample, lcd = mle_LogConcDEAD,
                           uselog = TRUE))

  density_df[sim == sim_val,
             mean_loglik_LogConcDEAD := mean_loglik_LogConcDEAD_eval]

}

# Save simulation density_df
data.table::fwrite(density_df, file = "sim_data/fig09_densities.csv")

# Save simulation points
data.table::fwrite(density_points, file = "sim_data/fig09_points.csv")
