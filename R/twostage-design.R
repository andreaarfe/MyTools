#' Evaluate operating characteristics of a two-stage design
#'
#' Computes type I error, power, and expected sample sizes under the null and
#' alternative hypotheses for a candidate two-stage single-arm design with a
#' binary endpoint.
#'
#' Decision rules:
#' - **Stage 1** (X1 ~ Binomial(n1, p)):
#'   - X1 <= r1: stop for futility (`r1 = -1` disables futility stopping).
#'   - X1 >= e1: declare efficacy (`e1 = n1 + 1` disables interim efficacy stop).
#'   - Otherwise: continue to stage 2.
#' - **Stage 2** (X2 ~ Binomial(n2, p), independent of X1):
#'   - X1 + X2 >= r: declare success.
#'
#' The irrevocability constraint `e1 >= r` guarantees that interim efficacy
#' cannot be overturned by stage 2 (worst case X2 = 0 still gives
#' X1 + X2 >= e1 >= r).
#'
#' @param n1 Integer. Stage-1 sample size.
#' @param n  Integer. Total sample size (n1 + n2).
#' @param r1 Integer. Futility boundary; stop if X1 <= r1. Use -1 for no
#'   futility stopping.
#' @param e1 Integer. Interim efficacy boundary; declare success if X1 >= e1.
#'   Use n1 + 1 for no interim efficacy stopping.
#' @param r  Integer. Final success threshold; success if X1 + X2 >= r.
#' @param p0 Numeric in (0, 1). Null response rate.
#' @param p1 Numeric in (0, 1). Alternative response rate (p1 > p0).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{n1, n, n2, r1, e1, r, p0, p1}{Design parameters (n2 = n - n1).}
#'     \item{alpha_actual}{Achieved type I error under p0.}
#'     \item{power_actual}{Achieved power under p1.}
#'     \item{en_null}{Expected sample size under p0.}
#'     \item{en_alt}{Expected sample size under p1.}
#'     \item{is_irrevocable}{Logical; TRUE when e1 >= r.}
#'   }
#' @examples
#' twostage_evaluate_design(n1 = 15L, n = 25L, r1 = 0L, e1 = 9L, r = 5L,
#'                          p0 = 0.10, p1 = 0.30)
#' @export
twostage_evaluate_design <- function(n1, n, r1, e1, r, p0, p1) {
  evaluate_design_cpp(as.integer(n1), as.integer(n), as.integer(r1),
                      as.integer(e1), as.integer(r),
                      as.numeric(p0), as.numeric(p1))
}
