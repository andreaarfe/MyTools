#' Filter designs to the Pareto frontier on sample size and expected sample sizes
#'
#' Retains only admissible (Pareto-optimal) designs with respect to three
#' criteria: maximum sample size (`n`), expected sample size under the null
#' (`en_null`), and expected sample size under the alternative (`en_alt`).
#' A design is admissible if no other design is weakly better on all three
#' axes with at least one strict improvement.
#'
#' @param designs A `data.frame` as returned by [find_feasible_designs()],
#'   containing at least columns `n`, `en_null`, and `en_alt`.
#'
#' @return A `data.frame` of admissible designs sorted by `en_null` ascending.
#'   Returns `designs` unchanged (empty) when the input has zero rows.
#' @export
find_admissible_designs <- function(designs) {
  if (nrow(designs) == 0L) return(designs)

  keep <- find_admissible_designs_cpp(as.integer(designs$n),
                                      as.numeric(designs$en_null),
                                      as.numeric(designs$en_alt))
  kept <- designs[keep, , drop = FALSE]
  kept[order(kept$en_null), , drop = FALSE]
}
