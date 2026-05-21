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
#'   overturned at the final analysis. Use `TRUE` when the interim efficacy
#'   stopping rule is **non-binding** (the trial may continue to stage 2 even
#'   after `X1 >= e1`), so that the final decision remains coherent with the
#'   interim declaration. Set to `FALSE` to include designs without this
#'   constraint.
#' @param simon Logical. If `TRUE`, restrict to Simon two-stage designs,
#'   i.e. designs with no interim stopping for efficacy (`e1 = n1 + 1`).
#'   Default `FALSE`. When `simon = TRUE` the `irrevocable` argument has no
#'   effect: because there is no interim efficacy declaration, the
#'   irrevocability constraint is vacuously satisfied regardless of `r`.
#'
#' @return A `data.frame` with one row per feasible design and columns
#'   `n1`, `n`, `n2`, `r1`, `e1`, `r`, `p0`, `p1`, `alpha_actual`,
#'   `power_actual`, `en_null`, `en_alt`, `is_irrevocable`.
#'   Returns an empty `data.frame` (with a message) when no designs are found.
#' @examples
#' find_feasible_designs(p0 = 0.10, p1 = 0.30, alpha = 0.05, power = 0.80,
#'                       n_max = 40L)
#' @export
find_feasible_designs <- function(p0, p1, alpha, power,
                                  n_max, n1_min = 5L,
                                  irrevocable = TRUE,
                                  simon = FALSE) {
  stopifnot(0 < p0, p0 < p1, p1 < 1,
            0 < alpha, alpha < 1,
            0 < power, power < 1,
            n_max >= n1_min + 1L)

  out <- find_feasible_designs_cpp(as.numeric(p0), as.numeric(p1),
                                   as.numeric(alpha), as.numeric(power),
                                   as.integer(n_max), as.integer(n1_min),
                                   as.logical(irrevocable),
                                   as.logical(simon))
  if (nrow(out) == 0L) {
    message("No feasible designs found. Try increasing n_max.")
    return(data.frame())
  }
  out
}
