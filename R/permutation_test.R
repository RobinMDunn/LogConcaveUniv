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
  dens_est <- LogConcDEAD::mlelcd(x = true_sample)

  # Sample n_obs from MLE log-concave density
  dens_sample <- LogConcDEAD::rlcd(n = n_obs, lcd = dens_est)

  # Combine all observations
  all_obs <- rbind(true_sample, dens_sample)

  # Initialize original sample test statistic
  orig_ts <- 0

  # Compute original sample test statistic
  for(i in 1:nrow(all_obs)) {

    # Get distance between all observations in true_sample and obs i
    P_n <- apply(matrix(true_sample, ncol = d), MARGIN = 1,
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

  # Reject H_0 if orig_ts > (B+1)*(1-alpha) quantile of shuffle_ts
  reject_null <- as.numeric(orig_ts > sort(shuffle_ts)[ceiling((B+1)*(1-alpha))])

  return(reject_null)
}
