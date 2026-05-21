#' Enumerate all feasible two-stage designs
#'
#' Exhaustively searches over all combinations of (n, n1, r, e1, r1) with
#' total sample size up to `n_max` and returns designs that satisfy the type I
#' error and power constraints.
#'
#' A design is feasible when:
#' - `alpha_actual <= alpha`
#' - `power_actual >= power`
#' - `r1 < e1` (non-trivial continuation region)
#' - `e1 >= r` (irrevocability, enforced only when `irrevocable = TRUE`)
#'
#' @param p0          Numeric in (0, 1). Null response rate.
#' @param p1          Numeric in (0, 1). Alternative response rate (p1 > p0).
#' @param alpha       Numeric in (0, 1). Maximum allowable type I error.
#' @param power       Numeric in (0, 1). Minimum required power.
#' @param n_max       Integer. Maximum total sample size to consider.
#' @param n1_min      Integer. Minimum stage-1 sample size (default 5).
#' @param irrevocable Logical. If `TRUE` (default), restrict to designs where
#'   `e1 >= r`, guaranteeing that an interim efficacy declaration cannot be
#'   overturned at the final analysis. Set to `FALSE` to include designs
#'   without this constraint.
#'
#' @return A `data.frame` with one row per feasible design and columns
#'   `n1`, `n`, `n2`, `r1`, `e1`, `r`, `p0`, `p1`, `alpha_actual`,
#'   `power_actual`, `en_null`, `en_alt`, `is_irrevocable`.
#'   Returns an empty `data.frame` (with a message) when no designs are found.
#' @export
find_feasible_designs <- function(p0, p1, alpha, power,
                                  n_max, n1_min = 5L,
                                  irrevocable = TRUE) {
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

        # Skip r = 0: every trial succeeds trivially (type I error = 1).
        # When irrevocable, also skip r > n1: e1 >= r would require e1 > n1,
        # which is impossible (e1 = n1+1 encodes no interim stop, not e1 > n1).
        if (r == 0L) next
        if (irrevocable && r > n1) next

        # e1 range: [r, n1+1] when irrevocable (enforces e1 >= r);
        #           [0, n1+1] otherwise.
        e1_min <- if (irrevocable) r else 0L
        for (e1 in e1_min:(n1 + 1L)) {

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
