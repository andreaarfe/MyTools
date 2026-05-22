#' Filter designs to the 3D Pareto frontier (internal)
#'
#' Retains only admissible (Pareto-optimal) designs with respect to three
#' criteria: maximum sample size (`n`), expected sample size under the null
#' (`en_null`), and expected sample size under the alternative (`en_alt`).
#' A design is admissible if no other design is weakly better on all three
#' axes with at least one strict improvement.
#'
#' Used internally by [twostage_admissible_designs()]. Not exported; users
#' should call [twostage_admissible_designs()] instead.
#'
#' @param designs A `data.frame` with at least columns `n`, `en_null`, and `en_alt`.
#' @return A `data.frame` of admissible designs sorted by `n` descending.
#'   Returns `designs` unchanged (empty) when the input has zero rows.
#' @keywords internal
twostage_find_admissible_designs <- function(designs) {
  if (nrow(designs) == 0L) return(designs)

  keep <- find_admissible_designs_cpp(as.integer(designs$n),
                                      as.numeric(designs$en_null),
                                      as.numeric(designs$en_alt))
  kept <- designs[keep, , drop = FALSE]
  rownames(kept) <- NULL
  kept[order(kept$n, decreasing = TRUE), , drop = FALSE]
}
