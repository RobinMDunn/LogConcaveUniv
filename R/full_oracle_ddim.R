#' Title
#'
#' @param data
#' @param B
#' @param alpha
#' @param mu
#' @param p
#' @param compute_ts
#'
#' @return
#' @export
full_oracle_ddim <- function(data, B, alpha, mu, p, compute_ts) {

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
      (1 - p) * mvtnorm::dmvnorm(x = Y_0, mean = rep(0, d), sigma = diag(d)) +
      p * mvtnorm::dmvnorm(x = Y_0, mean = 0 - mu, sigma = diag(d))

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
