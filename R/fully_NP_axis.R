#' Title
#'
#' @param data
#' @param B
#' @param alpha
#'
#' @return
#' @export
fully_NP_axis <- function(data, B, alpha) {

  # Extract number of observations and dimension
  n_obs <- nrow(data)
  d <- ncol(data)

  # Matrix of subsampled test statistics
  ts_mat <- matrix(NA, nrow = B, ncol = d)

  # Consider each dimension separately
  for(d_val in 1:d) {

    # Repeatedly subsample to get subsampled test statistic on each dimension
    for(b in 1:B) {

      # Split Y into Y_0 and Y_1
      Y_0_indices <- sample(1:n_obs, size = n_obs/2)

      Y_1_indices <- setdiff(1:n_obs, Y_0_indices)

      Y_0 <- matrix(data[Y_0_indices, ], ncol = d)

      Y_1 <- matrix(data[Y_1_indices, ], ncol = d)

      # Get KDE estimate for each dimension on D_1.
      # kde_D1 <- ks::kde(x = Y_1[, d_val], eval.points = Y_0[, d_val])
      kde_D1 <- ks::kde(x = Y_1[, d_val], eval.points = Y_0[, d_val],
                        h = ks::hpi(x = Y_1[, d_val]), binned = FALSE)

      # Evaluate KDE on D_0
      eval_kde_D0 <- kde_D1$estimate

      # Get log-concave MLE on D_0
      log_concave_D0 <- logcondens::logConDens(x = Y_0[, d_val])

      # Evaluate log-concave MLE on D_0
      eval_log_concave_D0 <- logcondens::evaluateLogConDens(
        xs = Y_0[, d_val], res = log_concave_D0, which = 2)[, 3]

      # Store dimension d, subsample b test stat
      ts_mat[b, d_val] <-
        exp(sum(log(eval_kde_D0)) - sum(log(eval_log_concave_D0)))

      # Break if you would reject based on current info
      if(sum(apply(ts_mat, MARGIN = 2, FUN = function(x) sum(x, na.rm = TRUE)) >=
             B * d /alpha) >= 1 & b < B) {
        ts_mat[(b+1):B, d_val] <- 0
        break
      }

    }

  }

  # Get average subsampled test stat for each dimension
  test_stats <- apply(ts_mat, MARGIN = 2, FUN = mean)

  # Reject H_0 if avg_ts >= 1/alpha
  return(list(reject_null = as.numeric(sum(test_stats >= d/alpha) >= 1)))

}
