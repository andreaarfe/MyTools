#' Operating characteristic curve for a two-stage design
#'
#' Evaluates the operating characteristics of a two-stage design at each
#' element of a user-supplied vector of response probabilities \code{p}.
#'
#' @param n1 Stage-1 sample size.
#' @param n  Total sample size (\code{n > n1}).
#' @param r1 Futility boundary: stop for futility if \eqn{X_1 \le r_1}
#'   (\code{-1} = no futility stop).
#' @param e1 Interim efficacy boundary: declare success if \eqn{X_1 \ge e_1}
#'   (\code{n1 + 1} = no interim efficacy stop).
#' @param r  Final success threshold: success if \eqn{X_1 + X_2 \ge r}.
#' @param p  Numeric vector of response probabilities in \eqn{(0, 1)}.
#'
#' @return A \code{data.frame} with one row per element of \code{p} and columns:
#'   \describe{
#'     \item{p}{Response probability.}
#'     \item{p_success}{Probability of trial success (early efficacy or
#'       \eqn{X_1 + X_2 \ge r}).}
#'     \item{p_futility}{Probability of early stopping for futility
#'       (\eqn{X_1 \le r_1}).}
#'     \item{p_efficacy1}{Probability of early stopping for efficacy at the
#'       interim analysis (\eqn{X_1 \ge e_1}).}
#'     \item{en}{Expected sample size.}
#'   }
#' @examples
#' twostage_evaluate_design_curve(n1 = 15L, n = 25L, r1 = 0L, e1 = 9L, r = 5L,
#'                                p = seq(0.05, 0.50, by = 0.05))
#' @export
twostage_evaluate_design_curve <- function(n1, n, r1, e1, r, p) {
  evaluate_design_curve_cpp(
    as.integer(n1), as.integer(n), as.integer(r1),
    as.integer(e1), as.integer(r), as.numeric(p)
  )
}
