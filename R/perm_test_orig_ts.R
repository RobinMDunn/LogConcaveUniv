#' Compute original permutation test statistic
#'
#' @description This function computes the original permutation test statistic,
#' keeping separate the original sample and data sampled from log-concave MLE.
#' Suppose Y contains the data from the original sample (`data`) and Y*
#' contains the data sampled from the log-concave MLE (`dens_sample`).
#' Suppose A_0 is the set of all balls centered at a point in union(Y, Y*),
#' P_n(A) is the proportion of observations in ball A out of all observations
#' in Y, and P*_n(A) is the proportion of observations in ball A out of all
#' observations in Y*. The test statistic equals
#' \eqn{sup_{A in A_0} |P_n(A) - P*_n(A)|}.
#'
#' @param data \eqn{n x d} data frame containing iid original observations.
#' One row per observation.
#' We wish to test whether the underlying density is log-concave.
#' @param dens_sample \eqn{n x d} data frame containing observations sampled
#' from the log-concave MLE
#' @param all_obs \eq{2n x d} data frame containing original observations and
#' data sampled from the log-concave MLE. Use the command
#' `rbind(data, dens_sample)` to construct this data frame.
#' @param d Dimension of the original sample
#'
#' @return Numeric value equal to original permutation test statistic
#' @export
perm_test_orig_ts <- function(data, dens_sample, all_obs, d) {

  # Initialize original sample test statistic
  orig_ts <- 0

  # Compute original sample test statistic
  for(i in 1:nrow(all_obs)) {

    # Get distance between all observations in data and obs i
    P_n <- apply(matrix(data, ncol = d), MARGIN = 1,
                 FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))

    # Get distance between all observations in dens_sample and obs i
    P_n_star <- apply(matrix(dens_sample, ncol = d), MARGIN = 1,
                      FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))

    # Get proportion of observations in sphere w/ each radius of size dist.
    # Save test stat that maximizes absolute difference.
    for(dist in unique(c(P_n, P_n_star))) {
      ts <- abs(mean(P_n <= dist) - mean(P_n_star <= dist))
      if(ts > orig_ts) {
        orig_ts <- ts
      }
    }

  }

  # Return original test stat
  return(orig_ts)

}
