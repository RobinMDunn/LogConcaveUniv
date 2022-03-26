#' Partial oracle, random projection test for log-concavity
#'
#' @description Run the partial oracle, random projection test of
#' "H_0: true density is log-concave" versus
#' "H_1: true density is not log-concave."
#' This function computes the numerator density by fitting a mixture of
#' two d-dimensional Normal distributions. This is a "partial oracle"
#' because this approach uses the fact that the density is a mixture of
#' Normals, but it estimates all parameters of the mixture.
#' This method computes test statistics by projecting the data onto `n_proj`
#' random directions. We construct `B` subsampled test statistics on each of the
#' `n_proj` directions. The final test statistic is the average of the
#' `B` x `n_proj` test statistics. This method rejects H_0 if the average
#' exceeds 1/alpha.
#'
#' @param data \eqn{n x d} data frame containing iid observations.
#' One row per observation.
#' We wish to test whether the underlying density is log-concave.
#' @param B Number of repeated subsamples on each random projection
#' @param n_proj Number of random projections
#' @param alpha Significance level
#' @param compute_ts Indicator for whether to compute test statistic.
#' Set `compute_ts = 0` to stop early if rejection is guaranteed after some
#' b < B subsamples.
#' Set `compute_ts = 1` to perform all B subsamples and compute the
#' test statistic.
#'
#' @return List containing `test_stat` and `reject_null`.
#' \itemize{
#'   \item `test_stat` --- If `compute_ts = 1`, this is the final test statistic,
#'   averaged over `B` x `n_proj` subsamples. If `compute_ts = 0`, this is NA.
#'   \item `reject_null` --- Indicator that equals 1 if we reject H_0 at level
#'   `alpha` and 0 if we do not reject H_0 at level `alpha`.
#' }
#' @export
partial_oracle_randproj <- function(data, B, n_proj, alpha, compute_ts) {

  # Extract number of observations and dimension
  n_obs <- nrow(data)
  d <- ncol(data)

  # Matrix of subsampled test statistics
  ts_mat <- matrix(NA, nrow = n_proj, ncol = B)

  # Repeatedly get test statistic on different projections
  for(proj in 1:n_proj) {

    # Get random vector for projection
    random_vector <- rnorm(n = d, mean = 0, sd = 1)

    random_vector <- random_vector / sqrt(sum(random_vector^2))

    for(b in 1:B) {

      # Split Y into Y_0 and Y_1
      Y_0_indices <- sample(1:n_obs, size = n_obs/2)

      Y_1_indices <- setdiff(1:n_obs, Y_0_indices)

      Y_0 <- matrix(data[Y_0_indices, ], ncol = d)

      Y_1 <- matrix(data[Y_1_indices, ], ncol = d)

      # Get projections of Y_0 and Y_1
      Y_0 <- as.numeric(Y_0 %*% random_vector)

      Y_1 <- as.numeric(Y_1 %*% random_vector)

      # Remove previously fitted mclust_dens_D1 object
      if(exists("mclust_dens_D1")) { rm(mclust_dens_D1) }

      # Get two-component Gaussian mixture on D_1.
      ptm <- proc.time()
      mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1, G = 2,
                                                  modelNames = "V",
                                                  warn = FALSE, verbose = FALSE,
                                                  plot = FALSE))
      proc.time() - ptm
      if(is(mclust_dens_D1, "try-error") | is.null(mclust_dens_D1)) {
        mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1, G = 2,
                                                    modelNames = "E",
                                                    warn = FALSE,
                                                    verbose = FALSE,
                                                    plot = FALSE))
      }

      # If mclust_dens_D1 is not an Mclust object, fit a single Normal density.
      # (This is very rare in simulations.)
      if(!("Mclust" %in% class(mclust_dens_D1))) {
        print("Fitting single Normal density")
        mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1, G = 1,
                                                    modelNames = "V",
                                                    warn = FALSE,
                                                    verbose = FALSE,
                                                    plot = FALSE))
      }

      # Evaluate Gaussian mixture on D_0
      eval_gauss_mix_D0 <- mclust::predict.densityMclust(object = mclust_dens_D1,
                                                         newdata = Y_0,
                                                         what = "dens")

      # Get log-concave MLE on D_0
      log_concave_D0 <- logcondens::logConDens(x = Y_0, smoothed = FALSE)

      # Evaluate log-concave MLE on D_0
      eval_log_concave_D0 <- logcondens::evaluateLogConDens(
        xs = Y_0, res = log_concave_D0, which = 2)[, 3]

      # Store test stat for the projection
      ts_mat[proj, b] <-
        exp(sum(log(eval_gauss_mix_D0)) - sum(log(eval_log_concave_D0)))

      # Stop if test stat is NA
      stopifnot(!is.na(ts_mat[proj, b]))

    }

    # Stop early if likelihood ratio already guarantees rejection
    if(compute_ts == 0 & sum(ts_mat, na.rm = TRUE) >= n_proj * B / alpha) {
      ts_mat[is.na(ts_mat)] <- 0
      test_stat <- NA_real_
      reject_null <- 1
      break
    }

  }

  # If not computing test stat but cannot reject, set test_stat to NA and
  # reject_null to 0
  if(compute_ts == 0 & sum(ts_mat) < n_proj * B / alpha) {
    test_stat <- NA_real_
    reject_null <- 0
  }

  # If computing test stat, get average subsampled test stat across projections
  if(compute_ts == 1) {
    test_stat <- mean(ts_mat)
    reject_null <- as.numeric(test_stat >= 1/alpha)
  }

  # Return test statistic and whether to reject H_0.
  # Reject H_0 if avg_ts >= 1/alpha
  return(list(test_stat = test_stat,
              reject_null = reject_null))

}
