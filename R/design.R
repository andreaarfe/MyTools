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
#' @export
evaluate_design <- function(n1, n, r1, e1, r, p0, p1) {
  n2 <- n - n1

  .oc <- function(p) {
    # Full PMF of X1 (indexed 1-based: pmf_x1[k+1] = P(X1 = k))
    pmf_x1 <- dbinom(0:n1, size = n1, prob = p)

    # P(interim efficacy): X1 in {e1, ..., n1}
    p_eff1 <- if (e1 > n1) 0 else sum(pmf_x1[(e1 + 1L):(n1 + 1L)])

    # P(futility stop): X1 in {0, ..., r1}
    p_fut1 <- if (r1 < 0L) 0 else sum(pmf_x1[1L:(r1 + 1L)])

    # Stage-2 success contribution: X1 in {r1+1, ..., e1-1}.
    # When r1 = e1 - 1 the continuation region is empty; guard explicitly
    # because seq.int(a, b) with a > b produces a *decreasing* sequence
    # rather than an empty one, which would double-count boundary values.
    cont_vals <- if (r1 + 1L <= e1 - 1L) {
      seq.int(r1 + 1L, e1 - 1L)
    } else {
      integer(0L)
    }
    p_success_2 <- if (length(cont_vals) == 0L) {
      0
    } else {
      sum(vapply(cont_vals, function(x1) {
        need <- r - x1
        p_x2 <- if (need <= 0L) {
          1
        } else if (need > n2) {
          0
        } else {
          pbinom(need - 1L, size = n2, prob = p, lower.tail = FALSE)
        }
        pmf_x1[x1 + 1L] * p_x2
      }, numeric(1L)))
    }

    p_success      <- p_eff1 + p_success_2
    p_stop_interim <- p_fut1 + p_eff1
    en             <- n1 * p_stop_interim + n * (1 - p_stop_interim)
    list(p_success = p_success, en = en)
  }

  oc0 <- .oc(p0)
  oc1 <- .oc(p1)

  list(
    n1             = n1,
    n              = n,
    n2             = n2,
    r1             = r1,
    e1             = e1,
    r              = r,
    p0             = p0,
    p1             = p1,
    alpha_actual   = oc0$p_success,
    power_actual   = oc1$p_success,
    en_null        = oc0$en,
    en_alt         = oc1$en,
    is_irrevocable = (e1 >= r)
  )
}
