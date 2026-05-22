#' Find all admissible two-stage designs
#'
#' Combines [twostage_find_feasible_designs()] and
#' [twostage_find_admissible_designs()] into a single call. Returns the
#' Pareto-optimal designs on three criteria —
#' maximum sample size (`n`), expected sample size under the null (`en_null`),
#' and expected sample size under the alternative (`en_alt`) — and labels the
#' **minimax** (minimum `n`) and **optimal** (minimum `en_null`) designs.
#'
#' @param p0          Numeric in (0, 1). Null response rate.
#' @param p1          Numeric in (0, 1). Alternative response rate (p1 > p0).
#' @param alpha       Numeric in (0, 1). Maximum allowable type I error.
#' @param power       Numeric in (0, 1). Minimum required power.
#' @param n_max       Integer. Maximum total sample size to consider (default 50).
#' @param n1_min      Integer. Minimum stage-1 sample size (default 1).
#' @param irrevocable Logical. If `TRUE`, restrict to designs where `e1 >= r`,
#'   so that an interim efficacy declaration cannot be overturned at the final
#'   analysis. Use `TRUE` when the interim efficacy stopping rule is
#'   **non-binding** (the trial may continue to stage 2 even after `X1 >= e1`),
#'   so that the final decision remains coherent with the interim declaration.
#'   Default `FALSE`.
#' @param simon Logical. If `TRUE`, restrict to Simon two-stage designs,
#'   i.e. designs with no interim stopping for efficacy (`e1 = n1 + 1`).
#'   Default `FALSE`. When `simon = TRUE` the `irrevocable` argument has no
#'   effect: because there is no interim efficacy declaration, the
#'   irrevocability constraint is vacuously satisfied regardless of `r`.
#'
#' @return A `data.frame` with one row per admissible design and columns
#'   `n`, `n1`, `n2`, `r1`, `e1`, `r`, `alpha_actual`,
#'   `power_actual`, `en_null`, `en_alt`, `is_irrevocable`, and `design_type`.
#'   The `design_type` column is a character vector with values:
#'   \describe{
#'     \item{`"minimax"`}{The admissible design with the smallest `n`
#'       (ties broken by `en_null`).}
#'     \item{`"optimal"`}{The admissible design with the smallest `en_null`
#'       (ties broken by `n`).}
#'     \item{`"minimax, optimal"`}{A design that is simultaneously minimax and
#'       optimal.}
#'     \item{`""`}{All other admissible designs.}
#'   }
#'   Returns an empty `data.frame` when no feasible designs are found.
#'
#' @examples
#' twostage_admissible_designs(p0 = 0.10, p1 = 0.30, alpha = 0.05, power = 0.80,
#'                             n_max = 40L)
#' @export
twostage_admissible_designs <- function(p0, p1, alpha, power,
                                        n_max = 50L, n1_min = 1L,
                                        irrevocable = FALSE,
                                        simon = FALSE) {
  feasible  <- twostage_find_feasible_designs(p0, p1, alpha, power,
                                              n_max, n1_min, irrevocable, simon)
  if (nrow(feasible) == 0L) return(feasible)

  admissible <- twostage_find_admissible_designs(feasible)
  if (nrow(admissible) == 0L) return(admissible)

  # Minimax: smallest n, tie-break by en_null
  mm_idx <- which(admissible$n == min(admissible$n))
  mm_idx <- mm_idx[which.min(admissible$en_null[mm_idx])]

  # Optimal: smallest en_null, tie-break by n
  opt_idx <- which(admissible$en_null == min(admissible$en_null))
  opt_idx <- opt_idx[which.min(admissible$n[opt_idx])]

  design_type <- rep("", nrow(admissible))
  design_type[mm_idx]  <- "minimax"
  design_type[opt_idx] <- "optimal"
  if (mm_idx == opt_idx) design_type[mm_idx] <- "minimax, optimal"

  admissible$design_type <- design_type
  admissible
}
