# Thin named wrappers around base R's dbinom/pbinom to make intent explicit.

# P(X = k), X ~ Binom(n, p)
binom_pmf <- function(k, n, p) dbinom(k, size = n, prob = p)

# P(X <= k), X ~ Binom(n, p)
binom_cdf <- function(k, n, p) pbinom(k, size = n, prob = p)

# P(X >= k) = 1 - P(X <= k-1), X ~ Binom(n, p)
binom_upper <- function(k, n, p) {
  pbinom(k - 1L, size = n, prob = p, lower.tail = FALSE)
}
