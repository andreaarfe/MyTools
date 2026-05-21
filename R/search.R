# Enumerate all feasible two-stage designs with total sample size n <= n_max.
#
# A design is feasible when:
#   - alpha_actual <= alpha
#   - power_actual >= power
#   - e1 >= r         (irrevocability constraint)
#   - r1 < e1         (non-trivial continuation region)
#
# Parameters:
#   p0, p1   : null and alternative response rates (0 < p0 < p1 < 1)
#   alpha    : maximum allowable type I error
#   power    : minimum required power
#   n_max    : maximum total sample size
#   n1_min   : minimum stage-1 sample size (default 5)
#
# Returns a data.frame with one row per feasible design and columns:
#   n1, n, n2, r1, e1, r, p0, p1,
#   alpha_actual, power_actual, en_null, en_alt, is_irrevocable
find_feasible_designs <- function(p0, p1, alpha, power,
                                  n_max, n1_min = 5L) {
  stopifnot(0 < p0, p0 < p1, p1 < 1,
            0 < alpha, alpha < 1,
            0 < power, power < 1,
            n_max >= n1_min + 1L)

  results <- vector("list", 1024L)  # pre-allocate; will grow if needed
  count   <- 0L

  for (n in (n1_min + 1L):n_max) {
    for (n1 in n1_min:(n - 1L)) {
      n2 <- n - n1

      # r ranges over possible final success thresholds
      for (r in 0L:n) {

        # Irrevocability: e1 >= r, but e1 <= n1, so impossible when r > n1.
        # Also skip r = 0: every trial succeeds trivially (type I error = 1).
        if (r == 0L || r > n1) next

        # e1 in [r, n1+1]; n1+1 encodes "no interim efficacy stopping"
        for (e1 in r:(n1 + 1L)) {

          # r1 in {-1, 0, ..., e1-1}; -1 encodes "no futility stopping"
          for (r1 in (-1L):(e1 - 1L)) {

            d <- evaluate_design(n1, n, r1, e1, r, p0, p1)

            if (d$alpha_actual <= alpha && d$power_actual >= power) {
              count <- count + 1L
              if (count > length(results))
                length(results) <- 2L * length(results)
              results[[count]] <- d
            }
          }
        }
      }
    }
  }

  if (count == 0L) {
    message("No feasible designs found. Try increasing n_max.")
    return(data.frame())
  }

  do.call(rbind, lapply(results[seq_len(count)], as.data.frame))
}
