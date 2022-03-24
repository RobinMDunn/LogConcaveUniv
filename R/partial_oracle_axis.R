#' Partial oracle, axis-aligned test for log-concavity
#'
#' @description Run the partial oracle axis-aligned test of
#' "H_0: true density is log-concave" versus
#' "H_1: true density is not log-concave."
#' This function computes the numerator density by fitting a mixture of
#' two d-dimensional Normal distributions. This is a "partial oracle"
#' because this approach uses the fact that the density is a mixture of
#' Normals, but it estimates all parameters of the mixture.
#' This method computes test statistics (averaged over B subsamples) on each of
#' the d axis directions. The test rejects H_0 if at least one of the test
#' statistics exceeds d/alpha.
#'
#' @param data \eqn{n x d} data frame containing iid observations.
#' One row per observation.
#' We wish to test whether the underlying density is log-concave.
#' @param B Number of repeated subsamples for test statistic construction
#' @param alpha Significance level
#'
#' @return List containing `reject_null`.
#' \itemize{
#'   \item `reject_null` --- Indicator that equals 1 if we reject H_0 at level
#'   `alpha` and 0 if we do not reject H_0 at level `alpha`.
#' }
#' @export
partial_oracle_axis <- function(data, B, alpha) {

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

      # Remove previously fitted mclust_dens_D1 object
      if(exists("mclust_dens_D1")) { rm(mclust_dens_D1) }

      # Get two-component Gaussian mixture estimate for each dimension on D_1.
      # Try unequal variance model. If an error, use equal variance model.
      # (Any density is valid.)
      # modelNames:
      #   - V: unequal variance
      #   - E: equal variance

      mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1[, d_val], G = 2,
                                                  modelNames = "V",
                                                  warn = FALSE, verbose = FALSE,
                                                  plot = FALSE))
      if(is(mclust_dens_D1, "try-error") | is.null(mclust_dens_D1)) {
        mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1[, d_val], G = 2,
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
                                                         newdata = Y_0[, d_val],
                                                         what = "dens")


      # Stop if max evaluated density is over 100 (probably convergence issue)
      # (Never occurs in simulations.)
      stopifnot(max(eval_gauss_mix_D0) < 100)

      # Get log-concave MLE on D_0
      log_concave_D0 <- logcondens::logConDens(x = Y_0[, d_val])

      # Evaluate log-concave MLE on D_0
      eval_log_concave_D0 <- logcondens::evaluateLogConDens(
        xs = Y_0[, d_val], res = log_concave_D0, which = 2)[, 3]

      # Store dimension d, subsample B test stat
      ts_mat[b, d_val] <-
        exp(sum(log(eval_gauss_mix_D0)) - sum(log(eval_log_concave_D0)))

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

  # Reject H_0 if some dimension's test_stat >= d/alpha
  return(list(reject_null = as.numeric(sum(test_stats >= d/alpha) >= 1)))

}
