#' Binomial probability mass function
#'
#' @param k Integer. Number of successes.
#' @param n Integer. Number of trials.
#' @param p Numeric in (0, 1). Success probability.
#' @return P(X = k) where X ~ Binomial(n, p).
#' @keywords internal
binom_pmf <- function(k, n, p) dbinom(k, size = n, prob = p)

#' Binomial cumulative distribution function
#'
#' @param k Integer. Number of successes.
#' @param n Integer. Number of trials.
#' @param p Numeric in (0, 1). Success probability.
#' @return P(X <= k) where X ~ Binomial(n, p).
#' @keywords internal
binom_cdf <- function(k, n, p) pbinom(k, size = n, prob = p)

#' Binomial upper tail probability
#'
#' Computes P(X >= k) using `pbinom(..., lower.tail = FALSE)` for numerical
#' stability.
#'
#' @param k Integer. Lower bound (inclusive) for the upper tail.
#' @param n Integer. Number of trials.
#' @param p Numeric in (0, 1). Success probability.
#' @return P(X >= k) where X ~ Binomial(n, p).
#' @keywords internal
binom_upper <- function(k, n, p) {
  pbinom(k - 1L, size = n, prob = p, lower.tail = FALSE)
}
