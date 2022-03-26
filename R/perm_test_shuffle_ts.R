#' Compute shuffled permutation test statistics
#'
#' @description This function computes B shuffled permutation test statistics,
#' with observations shuffled between the original sample and the synthetic
#' data, sampled from the log-concave MLE. For each test statistic, suppose
#' the observed and synthetic data are randomly partitioned into Y(b) and
#' Y*(b).
#' Suppose A_0 is the set of all balls centered at a point in union(Y(b), Y*(b)),
#' \eqn{P_{n,b}(A)} is the proportion of observations in ball A out of all
#' observations in Y(b), and \eqn{P*_{n,b}(A)} is the proportion of observations in
#' ball A out of all observations in Y*(b). The shuffled test statistic equals
#' \eqn{sup_{A in A_0} |P_{n,b}(A) - P*_{n,b}(A)|}.
#'
#' @param all_obs \eqn{2n x d} data frame containing original observations and
#' data sampled from the log-concave MLE
#' @param B Number of shuffled test statistics to compute
#' @param d Dimension of an observation in the original sample
#' @param n_obs Number of observations in the original sample
#'
#' @return Numeric vector containing shuffled permutation test statistics
#' @export
perm_test_shuffle_ts <- function(all_obs, B, d, n_obs) {

  # Initialize permutation test statistics
  shuffle_ts <- rep(0, B)

  # Run permutation test
  for(b in 1:B) {

    # Shuffle order of observations
    index_shuffle <- sample(1:nrow(all_obs), size = nrow(all_obs),
                            replace = FALSE)

    # Compute test statistic for permutation
    for(i in 1:nrow(all_obs)) {

      # Get distance between all observations in first sample and obs i
      P_n <- apply(matrix(all_obs[index_shuffle[1:n_obs], ], ncol = d),
                   MARGIN = 1,
                   FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))

      # Get distance between all observations in second sample and obs i
      P_n_star <- apply(matrix(all_obs[index_shuffle[(n_obs+1):(2*n_obs)], ],
                               ncol = d),
                        MARGIN = 1,
                        FUN = function(x) sqrt(sum((all_obs[i, ] - x)^2)))

      # Get proportion of observations in sphere w/ each radius of size dist.
      # Save test stat that maximizes absolute difference.
      for(dist in unique(c(P_n, P_n_star))) {
        ts <- abs(mean(P_n <= dist) - mean(P_n_star <= dist))
        if(ts > shuffle_ts[b]) {
          shuffle_ts[b] <- ts
        }
      }

    }

  }

  # Return vector of shuffled test stats
  return(shuffle_ts)
}
