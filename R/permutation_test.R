#' Title
#'
#' @param data
#' @param B
#' @param alpha
#'
#' @return
#' @export
permutation_test <- function(data, B, alpha) {

  # Get MLE log-concave density
  dens_est <- LogConcDEAD::mlelcd(x = data)

  # Number of observations and dimension of original sample
  n_obs <- nrow(data)
  d <- ncol(data)

  # Sample n_obs from MLE log-concave density
  dens_sample <- LogConcDEAD::rlcd(n = n_obs, lcd = dens_est)

  # Combine all observations
  all_obs <- rbind(data, dens_sample)

  # Get original sample test statistic
  orig_ts <- perm_test_orig_ts(data = data, all_obs = all_obs, d = d,
                               dens_sample = dens_sample)

  # Get permutation test shuffled statistics
  shuffle_ts <- perm_test_shuffle_ts(all_obs = all_obs, B = B, d = d,
                                     n_obs = n_obs)

  # Reject H_0 if orig_ts > (B+1)*(1-alpha) quantile of shuffle_ts
  reject_null <- as.numeric(orig_ts > sort(shuffle_ts)[ceiling((B+1)*(1-alpha))])

  # Return orig test stat, vector of shuffled test stats, rejection decision
  return(list(orig_ts = orig_ts,
              shuffle_ts = shuffle_ts,
              reject_null = reject_null))
}
