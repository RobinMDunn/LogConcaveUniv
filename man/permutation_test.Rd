% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutation_test.R
\name{permutation_test}
\alias{permutation_test}
\title{Permutation test for log-concavity}
\usage{
permutation_test(data, B, alpha)
}
\arguments{
\item{data}{\eqn{n x d} data frame containing iid observations.
One row per observation.
We wish to test whether the underlying density is log-concave.}

\item{B}{Number of shuffled test statistics to compute}

\item{alpha}{Significance level}
}
\value{
List containing \code{orig_ts}, \code{shuffle_ts}, and \code{reject_null}.
\itemize{
\item \code{orig_ts} --- Test statistic on sample that keeps separate the
original sample and data sampled from log-concave MLE.
\item \code{shuffle_ts} --- Vector of B test statistics, computed on shufflings
of the original sample and data sampled from log-concave MLE.
\item \code{reject_null} --- Indicator that equals 1 if we reject H_0 at level
\code{alpha} and 0 if we do not reject H_0 at level \code{alpha}.
}
}
\description{
Run the permutation test of
"H_0: true density is log-concave" versus
"H_1: true density is not log-concave."
This method computes the log-concave MLE on the original sample,
draws a synthetic sample from the log-concave MLE, and computes an original
test statistic based on the difference (between the original sample and
the synthetic sample) in proportions of observations contained in certain
spheres. This method also computes B statistics with observations shuffled
between the original sample and the synthetic sample. The test rejects H_0
if the original test statistic exceeds the (B+1)*(1-alpha) quantile of
the shuffled test statistics.
}
