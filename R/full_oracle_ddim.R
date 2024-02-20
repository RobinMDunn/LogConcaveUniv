#' Full oracle d-dimensional test for log-concavity
#'
#' @description Run the full oracle d-dimensional test of
#' "H_0: true density is log-concave" versus
#' "H_1: true density is not log-concave."
#' This function assumes that the underlying distribution is a d-dimensional
#' mixture of two Normal distributions with the form
#' \eqn{(1-p) N(0, I_d) + p N(-mu, sigma^2 I_d)}. In addition, this approach assumes
#' that we know
#' the true underlying density. (Hence, this is a helpful theoretical
#' comparison, but it likely will not be used in practice.)
#' This method averages test statistics over B subsamples and
#' rejects H_0 if the average exceeds 1/alpha.
#'
#' @param data \eqn{n x d} data frame containing iid observations.
#' One row per observation.
#' We wish to test whether the underlying density is log-concave.
#' @param B Number of repeated subsamples for test statistic construction
#' @param alpha Significance level
#' @param mu Mean parameter from the underlying distribution
#' \eqn{(1-p) N(0, I_d) + p N(-mu, sigma^2 I_d)}
#' @param sigma Standard deviation parameter from the underlying distribution
#' \eqn{(1-p) N(0, I_d) + p N(-mu, sigma^2 I_d)}
#' @param p Mixing parameter from the underlying distribution
#' \eqn{(1-p) N(0, I_d) + p N(-mu, sigma^2 I_d)}
#' @param compute_ts Indicator for whether to compute test statistic.
#' Set `compute_ts = 0` to stop early if rejection is guaranteed after some
#' b < B subsamples.
#' Set `compute_ts = 1` to perform all B subsamples and compute the
#' test statistic.
#' @param normalizing_constant 1 divided by probability of region of truncation,
#' if distribution is truncated.
#'
#' @return List containing `test_stat` and `reject_null`.
#' \itemize{
#'   \item `test_stat` --- If `compute_ts = 1`, this is the final test statistic,
#'   averaged over B subsamples. If `compute_ts = 0`, this is NA.
#'   \item `reject_null` --- Indicator that equals 1 if we reject H_0 at level
#'   `alpha` and 0 if we do not reject H_0 at level `alpha`.
#' }
#' @export
full_oracle_ddim <- function(data, B, alpha, mu, sigma = 1, p, compute_ts,
                             normalizing_constant = 1) {

  # Extract number of observations and dimension
  n_obs <- nrow(data)
  d <- ncol(data)

  # Vector of subsampled test statistics
  ts_vec <- rep(NA, B)

  # Repeatedly subsample to get test statistic
  for(b in 1:B) {

    # Split Y into Y_0 and Y_1
    Y_0_indices <- sample(1:n_obs, size = floor(p * n_obs))

    Y_1_indices <- setdiff(1:n_obs, Y_0_indices)

    Y_0 <- matrix(data[Y_0_indices, ], ncol = d)

    Y_1 <- matrix(data[Y_1_indices, ], ncol = d)

    # Evaluate true density on D_0
    eval_true_D0 <-
      normalizing_constant *
      ((1 - p) * mvtnorm::dmvnorm(x = Y_0, mean = rep(0, d), sigma = diag(d)) +
       p * mvtnorm::dmvnorm(x = Y_0, mean = 0 - mu, sigma = sigma^2 * diag(d)))

    # Get log-concave MLE on D_0
    log_concave_D0 <- LogConcDEAD::mlelcd(x = Y_0)

    # Evaluate log-concave MLE on D_0
    eval_log_concave_D0 <- LogConcDEAD::dlcd(x = Y_0, lcd = log_concave_D0)

    # Compute test statistic
    ts_vec[b] <- exp(sum(log(eval_true_D0)) - sum(log(eval_log_concave_D0)))

    # Stop early if likelihood ratio already guarantees rejection
    if(compute_ts == 0 & sum(ts_vec[1:b]) >= B / alpha) {
      ts_vec[(b + 1):B] <- 0
      avg_ts <- NA_real_
      reject_null <- 1
      break
    }
  }

  # If not computing test stat but cannot reject, set test_stat to NA and
  # reject_null to 0
  if(compute_ts == 0 & sum(ts_vec) < B / alpha) {
    avg_ts <- NA_real_
    reject_null <- 0
  }

  # If computing test stat, get average subsampled test stat across projections.
  # Reject H_0 if avg_ts >= 1/alpha.
  if(compute_ts == 1) {
    avg_ts <- mean(ts_vec)
    reject_null <- as.numeric(avg_ts >= 1/alpha)
  }

  # Return test statistic and whether to reject H_0.
  return(list(test_stat = avg_ts,
              reject_null = reject_null))
}
