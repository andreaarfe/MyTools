#' Enumerate all feasible two-stage designs (internal)
#'
#' Exhaustively searches over all combinations of `(n, n1, r, e1, r1)` with
#' total sample size up to `n_max` and returns designs satisfying
#' `alpha_actual <= alpha`, `power_actual >= power`, `r1 < e1`, and (optionally)
#' the irrevocability constraint `e1 >= r`.
#'
#' Used internally by [twostage_admissible_designs()]. Not exported; users
#' should call [twostage_admissible_designs()] instead.
#'
#' @param p0,p1       Null and alternative response rates in `(0, 1)`, `p0 < p1`.
#' @param alpha,power Maximum type I error and minimum power, in `(0, 1)`.
#' @param n_max       Integer. Maximum total sample size to consider (default 50).
#' @param n1_min      Integer. Minimum stage-1 sample size (default 1).
#' @param irrevocable Logical. If `TRUE`, restrict to designs with `e1 >= r`.
#' @param simon       Logical. If `TRUE`, restrict to Simon designs (`e1 = n1 + 1`,
#'   no interim efficacy stop).
#' @return A `data.frame` with one row per feasible design and columns
#'   `n`, `n1`, `n2`, `r1`, `e1`, `r`, `alpha_actual`, `power_actual`,
#'   `en_null`, `en_alt`, `is_irrevocable`. Returns an empty `data.frame`
#'   (with a message) when no designs are found.
#' @keywords internal
twostage_find_feasible_designs <- function(p0, p1, alpha, power,
                                           n_max = 50L, n1_min = 1L,
                                           irrevocable = FALSE,
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
  out <- out[, c("n", "n1", "n2", "r1", "e1", "r",
                 "alpha_actual", "power_actual", "en_null", "en_alt",
                 "is_irrevocable"), drop = FALSE]
  rownames(out) <- NULL
  out
}
