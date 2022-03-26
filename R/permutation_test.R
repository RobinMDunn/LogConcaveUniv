#' Permutation test for log-concavity
#'
#' @description Run the permutation test of
#' "H_0: true density is log-concave" versus
#' "H_1: true density is not log-concave."
#' This method computes the log-concave MLE on the original sample,
#' draws a synthetic sample from the log-concave MLE, and computes an original
#' test statistic based on the difference (between the original sample and
#' the synthetic sample) in proportions of observations contained in certain
#' spheres. This method also computes B statistics with observations shuffled
#' between the original sample and the synthetic sample. The test rejects H_0
#' if the original test statistic exceeds the (B+1)*(1-alpha) quantile of
#' the shuffled test statistics.
#'
#' @param data \eqn{n x d} data frame containing iid observations.
#' One row per observation.
#' We wish to test whether the underlying density is log-concave.
#' @param B Number of shuffled test statistics to compute
#' @param alpha Significance level
#'
#' @return List containing `orig_ts`, `shuffle_ts`, and `reject_null`.
#' \itemize{
#'   \item `orig_ts` --- Test statistic on sample that keeps separate the
#'   original sample and data sampled from log-concave MLE.
#'   \item `shuffle_ts` --- Vector of B test statistics, computed on shufflings
#'   of the original sample and data sampled from log-concave MLE.
#'   \item `reject_null` --- Indicator that equals 1 if we reject H_0 at level
#'   `alpha` and 0 if we do not reject H_0 at level `alpha`.
#' }
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
