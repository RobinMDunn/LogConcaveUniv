#' Title
#'
#' @param all_obs
#' @param B
#' @param d
#' @param n_obs
#'
#' @return
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
