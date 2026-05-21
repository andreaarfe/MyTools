# Evaluate operating characteristics for a candidate two-stage design.
#
# Decision rules:
#   Stage 1 (observe X1 ~ Binom(n1, p)):
#     X1 <= r1  => stop for futility   (r1 = -1 means no futility stopping)
#     X1 >= e1  => declare efficacy    (e1 = n1+1 means no interim efficacy stop)
#     otherwise => continue to stage 2
#   Stage 2 (observe X2 ~ Binom(n2, p), independent of X1):
#     X1 + X2 >= r  => success
#
# Irrevocability: e1 >= r ensures that when X1 >= e1 and X2 = 0 (worst case),
# the total X1 + X2 >= e1 >= r, so interim efficacy cannot be overturned.
#
# Returns a named list with all design parameters and computed characteristics.
evaluate_design <- function(n1, n, r1, e1, r, p0, p1) {
  n2 <- n - n1

  .oc <- function(p) {
    # Full PMF of X1 (indexed 1-based: pmf_x1[k+1] = P(X1 = k))
    pmf_x1 <- dbinom(0:n1, size = n1, prob = p)

    # P(interim efficacy): X1 in {e1, ..., n1}
    p_eff1 <- if (e1 > n1) 0 else sum(pmf_x1[(e1 + 1L):(n1 + 1L)])

    # P(futility stop): X1 in {0, ..., r1}
    p_fut1 <- if (r1 < 0L) 0 else sum(pmf_x1[1L:(r1 + 1L)])

    # Stage-2 success contribution: X1 in {r1+1, ..., e1-1}
    cont_vals <- seq.int(r1 + 1L, e1 - 1L)
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
