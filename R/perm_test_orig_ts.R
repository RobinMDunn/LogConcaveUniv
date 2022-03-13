#' Title
#'
#' @param data
#' @param all_obs
#' @param d
#' @param dens_sample
#'
#' @return
#' @export
perm_test_orig_ts <- function(data, all_obs, d, dens_sample) {

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
