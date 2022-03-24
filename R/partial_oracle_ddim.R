#' Partial oracle d-dimensional test for log-concavity
#'
#' @description Run the partial oracle d-dimensional test of
#' "H_0: true density is log-concave" versus
#' "H_1: true density is not log-concave."
#' This function computes the numerator density by fitting a mixture of
#' two d-dimensional Normal distributions. This is a "partial oracle"
#' because this approach uses the fact that the density is a mixture of
#' Normals, but it estimates all parameters of the mixture.
#' This method averages test statistics over B subsamples and
#' rejects H_0 if the average exceeds 1/alpha.
#'
#' @param data \eqn{n x d} data frame containing iid observations.
#' One row per observation.
#' We wish to test whether the underlying density is log-concave.
#' @param B Number of repeated subsamples for test statistic construction
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
#'   averaged over B subsamples. If `compute_ts = 0`, this is NA.
#'   \item `reject_null` --- Indicator that equals 1 if we reject H_0 at level
#'   `alpha` and 0 if we do not reject H_0 at level `alpha`.
#' }
#' @export
partial_oracle_ddim <- function(data, B, alpha, compute_ts) {

  # Extract number of observations and dimension
  n_obs <- nrow(data)
  d <- ncol(data)

  # Vector of subsampled test statistics
  ts_vec <- rep(NA, B)

  # Repeatedly subsample to get test statistic
  for(b in 1:B) {

    # Split Y into Y_0 and Y_1
    Y_0_indices <- sample(1:n_obs, size = n_obs/2)

    Y_1_indices <- setdiff(1:n_obs, Y_0_indices)

    Y_0 <- matrix(data[Y_0_indices, ], ncol = d)

    Y_1 <- matrix(data[Y_1_indices, ], ncol = d)

    # Remove previously fitted mclust_dens_D1 object
    if(exists("mclust_dens_D1")) { rm(mclust_dens_D1) }

    # Get two-component Gaussian mixture estimate on D_1.
    # Try unequal variance model. If an error, use equal variance model.
    # (Any density is valid.)
    # modelNames:
    #   - V: unequal variance
    #   - E: equal variance
    #   - VVV: unstructured (other than symmetric) covariance matrices.
    #   - VEV: varying volume, equal shape, varying orientation (Banfield & Raftery, 1993)
    if(d == 1) {

      mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1, G = 2,
                                                  modelNames = "V",
                                                  warn = FALSE,
                                                  verbose = FALSE,
                                                  plot = FALSE))
      if(is(mclust_dens_D1, "try-error") | is.null(mclust_dens_D1)) {
        mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1, G = 2,
                                                    modelNames = "E",
                                                    warn = FALSE,
                                                    verbose = FALSE,
                                                    plot = FALSE))
      }

    } else {

      mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1, G = 2,
                                                  modelNames = "VVV",
                                                  warn = FALSE, verbose = FALSE,
                                                  plot = FALSE))
      if(is(mclust_dens_D1, "try-error") | is.null(mclust_dens_D1)) {
        mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1, G = 2,
                                                    modelNames = "VEV",
                                                    warn = FALSE,
                                                    verbose = FALSE,
                                                    plot = FALSE))
      }

    }

    # If mclust_dens_D1 is not an Mclust object, fit a single Normal density.
    # (This is very rare in simulations.)
    if(!("Mclust" %in% class(mclust_dens_D1))) {
      print("Fitting single Normal density")
      if(d == 1) {
        mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1, G = 1,
                                                    modelNames = "V",
                                                    warn = FALSE,
                                                    verbose = FALSE,
                                                    plot = FALSE))
      } else {
        mclust_dens_D1 <- try(mclust::densityMclust(data = Y_1, G = 1,
                                                    modelNames = "VVV",
                                                    warn = FALSE,
                                                    verbose = FALSE,
                                                    plot = FALSE))
      }
    }

    # # Parameters of Gaussian mixture
    # prop_est_D1 <- mclust_dens_D1$parameters$pro
    # mean_est_D1 <- mclust_dens_D1$parameters$mean (for 2 classes)
    # var_est_D1 <- mclust_dens_D1$parameters$variance$sigma # array of 2 d x d matrices for d > 1
    # var_est_D1 <- mclust_dens_D1$parameters$variance$sigmasq # vector length 2 for d = 1

    # Evaluate Gaussian mixture on D_0
    eval_gauss_mix_D0 <- mclust::predict.densityMclust(object = mclust_dens_D1,
                                                       newdata = Y_0,
                                                       what = "dens")

    # Stop if max evaluated density is over 100 (probably convergence issue)
    # (Never occurs in simulations.)
    stopifnot(max(eval_gauss_mix_D0) < 100)

    ## Same as
    # prop_est_D1[1]*mclust::dmvnorm(data = Y_0, mean = mean_est_D1[,1],
    #                                sigma = var_est_D1[,,1]) +
    #   prop_est_D1[2]*mclust::dmvnorm(data = Y_0, mean = mean_est_D1[,2],
    #                                  sigma = var_est_D1[,,2])

    # Get log-concave MLE on D_0
    log_concave_D0 <- LogConcDEAD::mlelcd(x = Y_0)

    # Evaluate log-concave MLE on D_0
    eval_log_concave_D0 <- LogConcDEAD::dlcd(x = Y_0, lcd = log_concave_D0)

    # Compute test statistic
    ts_vec[b] <- exp(sum(log(eval_gauss_mix_D0)) - sum(log(eval_log_concave_D0)))

    # Stop early if likelihood ratio already guarantees rejection
    if(compute_ts == 0 & sum(ts_vec[1:b]) >= B / alpha) {
      ts_vec[(b + 1):B] <- 0
      test_stat <- NA_real_
      reject_null <- 1
      break
    }
  }

  # If not computing test stat but cannot reject, set test_stat to NA and
  # reject_null to 0
  if(compute_ts == 0 & mean(ts_vec) < 1 / alpha) {
    test_stat <- NA_real_
    reject_null <- 0
  }

  # If computing test stat, get average subsampled test stat across projections
  if(compute_ts == 1) {
    test_stat <- mean(ts_vec)
    reject_null <- as.numeric(test_stat >= 1 / alpha)
  }

  # Return test statistic and whether to reject H_0.
  # Reject H_0 if avg_ts >= 1/alpha
  return(list(test_stat = test_stat,
              reject_null = reject_null))

}
